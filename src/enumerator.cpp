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

#include "enumerator.hpp"
#include "thread.hpp"

#include <numeric>
#include <shared_mutex>

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
  log::print(oss.str());
}

void enumerator::draw(const std::vector<graph::vertex> &clq) const {
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

enumerator::enumerator(graph G) : cache(1'000'000) {
  const std::size_t vertex_count = G.adjacency().size();

  auto begin = std::chrono::high_resolution_clock::now();

  // compute the neighbourhood degrees of all vertices for which they shall be
  // sorted based on their degeneracy
  std::pmr::monotonic_buffer_resource memory_pool;
  std::pmr::unordered_map<graph::vertex, std::size_t> neighs_degree(
      &memory_pool);
  neighs_degree.reserve(vertex_count);
  std::for_each(
      std::execution::par_unseq, G.A.cbegin(), G.A.cend(),
      [&neighs_degree, &G](const graph::adjacency_map::value_type &e) {
        neighs_degree[e.first] =
            std::accumulate(e.second.cbegin(), e.second.cend(), 0ul,
                            [&G](std::size_t total_degrees, graph::vertex u) {
                              return G.neighbours(u).size() + total_degrees;
                            });
      });

  // switch the vertex container from std::unordered_map to std::vector
  V.resize(vertex_count);
  std::transform(std::execution::par_unseq, G.A.begin(), G.A.end(), V.begin(),
                 [](graph::adjacency_map::value_type &e) {
                   return std::make_pair(e.first, std::move(e.second));
                 });

  // sort vertices based on the degree of their adjacency as it was proven
  // that it makes colouring vertices optimal (as close to brute-forced)
  std::sort(std::execution::par_unseq, V.begin(), V.end(),
            [&](const adjacency_vector::value_type &v,
                const adjacency_vector::value_type &u) {
              return (u.second.size() < v.second.size() ||
                      (u.second.size() == v.second.size() &&
                       neighs_degree.at(u.first) < neighs_degree.at(v.first)));
            });

  // allocate containers and reverse mapping between vertices and keys to
  // enumerate neighbours based on the order of vertices

  std::pmr::unordered_map<graph::vertex, key> mapping(&memory_pool);
  mapping.reserve(vertex_count);

  B.resize(vertex_count);
  std::for_each(
      std::execution::par_unseq, B.begin(), B.end(),
      [vertex_count](std::vector<bool> &mtx) { mtx.resize(vertex_count); });
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
    std::transform(std::execution::par_unseq, neighs.begin(), neighs.end(),
                   A[v].begin(), [&v_mtx, &mapping](graph::vertex u) {
                     key u_key = mapping.at(u);
                     v_mtx[u_key] = true;
                     return u_key;
                   });
    // then sort them too based on their induced indices for optimal colouring
    std::sort(std::execution::par_unseq, A[v].begin(), A[v].end());
  }

  auto end = std::chrono::high_resolution_clock::now();
  log::info("Vertices were enumerated in",
            log::time_diff(begin, end, log::bold));
}

inline void enumerator::cache_hit_progress() const {
  static const std::size_t cache_portion = cache.capacity() * .05;
  std::size_t hits = cache_hits;
  if (hits >= cache_portion && hits % cache_portion == 0) {
    log::info(COL_GREEN, "cache total_hits", hits);
  }
}

// TODO refactor the colouring part
enumerator::sorted_keys
enumerator::greedy_colour_sort(std::vector<key> &&vertices) const {
  assert(not vertices.empty());

  std::vector<colour> colours;
  std::unique_lock<std::shared_mutex> cache_write_lock(cache_mtx,
                                                       std::defer_lock);
  thread::scope_dtor unlock_cache([&cache_write_lock]() {
    if (cache_write_lock.owns_lock()) {
      cache_write_lock.unlock();
    }
  });

  {
    std::shared_lock<std::shared_mutex> cache_read_lock(cache_mtx);
    if (auto precomputed = cache.get(vertices); precomputed.has_value()) {
      cache_hits++;
      cache_hit_progress();
      return sorted_keys(precomputed.value());
    }
  }

  cache_write_lock.lock();
  if (auto precomputed = cache.get(vertices); precomputed.has_value()) {
    cache_hits++;
    cache_hit_progress();
    return sorted_keys(precomputed.value());
  }

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

  // sort the vertices in-place (overriding the incoming std::vector)
  colours.reserve(vertices.size());
  vertices.clear(); // clearing up for in-place sorting (after reserving memory)
  for (colour col = 0; col < col_class.size(); ++col) {
    std::fill_n(std::back_inserter(colours), col_class[col].size(), col + 1);
    // moving/appending the vertices is cheaper than std::copy_n
    std::move(col_class[col].begin(), col_class[col].end(),
              std::back_inserter(vertices));
  }

  // FIXME greedy colouring sorting reorders the neighbours
  // FIXME neighbours reference is lost when moved from
  return sorted_keys(cache.set(vertices, colours));
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
