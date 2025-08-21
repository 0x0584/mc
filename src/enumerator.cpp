// enumerator.cpp
//
// Copyright (C) 2024  0x0584 (Anas)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
// USA.

#include "profiler.hpp"
#include <future>
#include <numeric>
#include <set>
#include <shared_mutex>

#include "enumerator.hpp"
#include "pq.hpp"

namespace mc {

void enumerator::print() const {
  std::ostringstream oss;
  oss << "Enumerated Vertices\n";
  for (const auto &[v, neighs] : V) {
    oss << v << " { ";
    for (graph::vertex u : neighs) {
      oss << u << " ";
    }
    oss << "}\n";
  }
  oss << "\n";
  logger::print(oss.str());
}

void enumerator::draw(
    const std::vector<graph::vertex> &clq) const { // FIXME: optimise this
  std::set<graph::vertex> clique(clq.begin(), clq.end());
  std::ofstream file((args::stdin ? "out" : args::filename) + ".dot");
  std::set<std::set<graph::vertex>> edges;
  std::mutex edges_mtx, file_mtx, clique_mtx;
  std::vector<std::future<void>> neighs_callbacks;
  neighs_callbacks.reserve(V.size());

  file << "digraph {\n"
          "ratio=fill; overlap=false;\n"
          "node [width=0.1 height=0.1 fontsize=8 shape=plain];\n"
          "edge [color=orange penwidth=0.1];\n";
  for (const auto &tmp : V) {
    {
      std::unique_lock lock(file_mtx);
      file << tmp.first << " [label=" << tmp.first << " ";
      if (std::unique_lock lock(clique_mtx); clique.count(tmp.first) > 0) {
        file << "shape=circle";
      }
      file << "];\n";
    }
    neighs_callbacks.emplace_back(std::async(
        std::launch::deferred,
        [&](graph::vertex v, const graph::neighbours_set &neighs) {
          std::ostringstream oss;
          for (const auto &u : neighs) {
            if (std::unique_lock lock(edges_mtx); edges.count({u, v}) > 0) {
              continue;
            } else {
              edges.emplace(std::set<graph::vertex>{u, v});
            }
            oss << v << " -> " << u << " [arrowhead=none ";
            if (std::unique_lock lock(clique_mtx);
                clique.count(u) > 0 && clique.count(v) > 0) {
              oss << " color=black penwidth=0.7";
            }
            oss << "];\n";
          }
          std::unique_lock lock(file_mtx);
          file << oss.str();
        },
        tmp.first, tmp.second));
  }
  std::for_each(neighs_callbacks.begin(), neighs_callbacks.end(),
                [](auto &fn) { fn.get(); });
  file << "}\n";
  file.close();
}

enumerator::enumerator(
    graph &G) { // TODO: merge enumeration with graph building
  const std::size_t vertex_count = G.adjacency().size();

  auto begin = std::chrono::high_resolution_clock::now();

  // compute the neighbourhood degrees of all vertices for which they shall be
  // sorted based on their degeneracy

  std::pmr::unordered_map<graph::vertex, std::size_t> neighs_degree;
  neighs_degree.reserve(vertex_count);
  std::for_each(
      G.A.cbegin(), G.A.cend(),
      [&neighs_degree, &G](const graph::adjacency_map::value_type &e) {
        neighs_degree[e.first] =
            std::accumulate(e.second.cbegin(), e.second.cend(), 0ul,
                            [&G](std::size_t total_degrees, graph::vertex u) {
                              return G.neighbours(u).size() + total_degrees;
                            });
      });

  logger::debug("Computed neighbourhood degrees");

  // switch the vertex container from std::unordered_map to std::vector
  V.resize(vertex_count);
  std::transform(G.A.begin(), G.A.end(), V.begin(),
                 [](graph::adjacency_map::value_type &e) {
                   return std::make_pair(e.first, std::move(e.second));
                 });

  // sort vertices based on the degree of their adjacency as it was proven
  // that it makes colouring vertices optimal (as close to brute-forced)
  std::sort(V.begin(), V.end(),
            [&](const adjacency_vector::value_type &v,
                const adjacency_vector::value_type &u) {
              return (u.second.size() < v.second.size() ||
                      (u.second.size() == v.second.size() &&
                       neighs_degree.at(u.first) < neighs_degree.at(v.first)));
            });

  logger::debug("Sorted neighbourhood degrees");

  // allocate containers and reverse mapping between vertices and keys to
  // enumerate neighbours based on the order of vertices

  std::pmr::unordered_map<graph::vertex, key> mapping;
  mapping.reserve(vertex_count);

  B.resize(vertex_count);
  std::for_each(B.begin(), B.end(), [vertex_count](std::vector<bool> &mtx) {
    mtx.resize(vertex_count);
  });
  A.resize(vertex_count);
  for (key v = 0; v < vertex_count; ++v) {
    const adjacency_vector::value_type &v_pair =
        V[v]; // a pair of vertex its neighbours
    mapping.emplace(v_pair.first, v);
    A[v].resize(v_pair.second.size());
  }

  // fill vertex adjacency and sort neighbours based on the induced order
  for (key v = 0; v < vertex_count; ++v) {
    const graph::neighbours_set &neighs = V[v].second;
    std::vector<bool> &v_mtx = B[v];
    // induce the indices of neighbours based on the mapping order
    std::transform(neighs.begin(), neighs.end(), A[v].begin(),
                   [&v_mtx, &mapping](graph::vertex u) {
                     key u_key = mapping.at(u);
                     v_mtx[u_key] = true;
                     return u_key;
                   });
    // then sort them too based on their induced indices for optimal colouring
    std::sort(A[v].begin(), A[v].end());
  }

  auto end = std::chrono::high_resolution_clock::now();
  logger::info("Vertices were enumerated in",
               logger::time_diff(begin, end, logger::bold));
}

enumerator::sorted_keys enumerator::greedy_colour_sort(
    std::vector<key> &&vertices,
    lru_cache<std::vector<enumerator::key>, std::vector<enumerator::colour>>
        &cache,
    std::size_t &cache_hits) const {
  static thread_local size_t cache_hits_mod = 10'000;
  if (auto precomputed = cache.get(vertices); precomputed.has_value()) {
    cache_hits++;
    if (cache_hits % cache_hits_mod == 0) {
      cache_hits_mod += .5 * cache_hits_mod;
      logger::warn("current cache hits", cache_hits);
    }
    return sorted_keys(precomputed.value());
  } else {
    auto result = greedy_colour_sort(std::move(vertices));
    return sorted_keys(cache.set(result.keys(), result.colours()));
  }
}

// TODO use dynamic colouring instead
enumerator::sorted_keys
enumerator::greedy_colour_sort(std::vector<key> &&vertices) const {
  assert(not vertices.empty());

  logger::warn("greedy start");
  auto start = std::chrono::system_clock::now();
  std::vector<colour> colours;

  // dividing vertices into colour classes based on their order of adjacency
  // appearance, then rearrange them based on colour priority
  //
  std::vector<std::vector<key>> col_class(vertices.size());
  //
  // since at most there will be memory allocations as many vertices, it is
  // fine to not reserve memory beforehand (as far as I had tested!)
  for (key v : vertices) {
    const std::vector<bool> &v_mtx = B.at(v);
    colour col = 0;
    //
    // since also, colour sorting is used to prune unnecessary branching when
    // seeking exactitude, thus it will practically decrease performance if
    // checking adjacent colours ran in parallel (again, as far as I had
    // tested)
    //
    // although, since using an adjacency matrix is optimal for probing, this
    // would be far less overhead expense that justifies the memory footprint
    while (std::any_of(col_class[col].begin(), col_class[col].end(),
                       [&v_mtx](key u) { return v_mtx.at(u); })) {
      col++;
    }
    col_class[col].emplace_back(v);
  }

  auto end = std::chrono::system_clock::now();
  logger::warn("greedy done in", logger::time_diff(start, end));

  // sort the vertices in-place (overriding the incoming std::vector)
  colours.reserve(vertices.size());
  vertices.clear(); // clearing up for in-place sorting (after reserving memory)
  for (colour col = 0; col < col_class.size(); ++col) {
    std::fill_n(std::back_inserter(colours), col_class[col].size(), col + 1);
    // moving/appending the vertices is cheaper than std::copy_n
    std::move(col_class[col].begin(), col_class[col].end(),
              std::back_inserter(vertices));
  }

  end = std::chrono::system_clock::now();
  logger::warn("greedy finished", logger::time_diff(start, end));

  return sorted_keys(std::make_pair(std::move(vertices), std::move(colours)));
}

enumerator::sorted_keys
enumerator::dsatur_colour_sort(std::vector<key> &&vertices) const {
  static gc gc;

  assert(not vertices.empty());

  profiler_start("mann45.prof");

  const auto num_vert = vertices.size();

  logger::debug("dsatur V=", num_vert);

  // --- NEW: Create a key-to-index mapping ---
  std::pmr::vector<key> index_to_key(num_vert, gc.get_allocator());
  std::pmr::unordered_map<key, std::size_t> key_to_index(gc.get_allocator());
  key_to_index.reserve(num_vert);

  for (std::size_t i = 0; i < num_vert; ++i) {
    key_to_index[vertices[i]] = i;
    index_to_key[i] = vertices[i];
  }
  // Note: The key_to_index map can be a temporary std::unordered_map
  // that is discarded after this loop if memory is a concern.
  // --- END NEW ---

  // --- NEW: Use vectors instead of unordered_map ---
  std::pmr::vector<colour> colours_arr(num_vert, gc.get_allocator());
  std::fill(colours_arr.begin(), colours_arr.end(),
            static_cast<colour>(-1)); // Sentinel for uncoloured
  std::pmr::vector<std::size_t> satur_arr(num_vert, gc.get_allocator());
  // --- END NEW ---

  auto sort_keys = [&](std::size_t a_idx, std::size_t b_idx) {
    return A[index_to_key[a_idx]].size() > A[index_to_key[b_idx]].size() &&
           satur_arr[a_idx] < satur_arr[b_idx];
  };

  utils::pq<std::size_t, decltype(sort_keys)> uncoloured(sort_keys);

  auto start = std::chrono::system_clock::now();
  for (std::size_t i = 0; i < num_vert; ++i) {
    uncoloured.emplace(i);
  }

  logger::warn("dsatur start");

  // The size of used_colours should be capped by the number of vertices,
  // as the maximum possible color is num_vert - 1.
  std::pmr::vector<bool> used_colours(num_vert, gc.get_allocator());

  while (!uncoloured.empty()) {
    std::size_t u_idx = uncoloured.top();
    key u = index_to_key[u_idx];
    uncoloured.pop();

    const auto &u_adj = A[u];
    std::fill(used_colours.begin(), used_colours.end(), false);

    for (key v : u_adj) {
      // Check if 'v' has been coloured
      if (key_to_index.count(v) &&
          colours_arr[key_to_index.at(v)] != static_cast<colour>(-1)) {
        used_colours[colours_arr[key_to_index.at(v)]] = true;
      }
    }

    colour col = 0;
    while (col < num_vert && used_colours[col]) {
      col++;
    }
    colours_arr[u_idx] = col;

    for (key v : u_adj) {
      std::size_t v_idx = key_to_index.at(v);
      if (colours_arr[v_idx] == static_cast<colour>(-1)) {
        // The vertex 'v' is not yet coloured, so its saturation isn't finalized
        continue;
      }

      std::size_t count = 0;
      for (key w : A[v]) {
        if (key_to_index.count(w) &&
            colours_arr[key_to_index.at(w)] != static_cast<colour>(-1)) {
          count++;
        }
      }
      satur_arr[v_idx] = count;
    }
  }

  auto end = std::chrono::system_clock::now();
  logger::warn("dsatur done in", logger::time_diff(start, end));

  std::vector<std::pair<key, colour>> keys_cols;
  for (std::size_t i = 0; i < num_vert; ++i) {
    keys_cols.emplace_back(index_to_key[i], colours_arr[i]);
  }

  std::sort(keys_cols.begin(), keys_cols.end(),
            [](const auto &a, const auto &b) { return a.second < b.second; });

  std::vector<key> keys;
  std::vector<colour> cols;

  for (auto &[k, col] : keys_cols) {
    keys.emplace_back(k);
    cols.emplace_back(col);
  }

  end = std::chrono::system_clock::now();
  logger::warn("dsatur finished", logger::time_diff(start, end));
  logger::debug("V=", keys.size());

  profiler_stop();
  return sorted_keys(std::make_pair(std::move(keys), std::move(cols)));
}

bool enumerator::is_clique(const std::vector<key> &clique) const {
  for (key v : clique) {
    for (key u : clique) {
      if (v != u) {
        if (const std::vector<key> &neighs = A.at(u);
            std::find(std::execution::par_unseq, neighs.cbegin(), neighs.cend(),
                      v) == neighs.cend()) {
          return false;
        }
      }
    }
  }
  return true;
}
} // namespace mc
