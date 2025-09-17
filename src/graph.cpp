// graph.cpp
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

#include "graph.hpp"
#include "core.hpp"

#include <algorithm>
#include <numeric>

#include "thread.hpp"
#include <execution>

namespace mc {
void graph::print() const {
  std::ostringstream oss;
  oss << " Graph (vertices=" << n_vertices << ", edges=" << n_edges << ")\n";
  for (key u = 0ul; u < n_vertices; ++u) {
    oss << to_vertex(u) << " { ";
    for (auto [it, last] = neighbours(u); it != last; it++) {
      auto v = *it;
      oss << to_vertex(v) << " ";
    }
    oss << "}\n";
  }
  logger::print(oss.str());
}

constexpr static unsigned MAX_DIGITS = 12;

static inline bool is_digit(char c) noexcept {
  return static_cast<unsigned>(c - '0') <= 9;
}

static inline unsigned digit(char c) noexcept {
  return static_cast<unsigned>(c - '0');
}

static inline bool is_delimiter(char c) noexcept {
  return c == feed::deli || c == feed::sep;
}

static inline bool
read_single_vertex(graph::vertex &w, feed::buffer_iterator &__restrict it,
                   const feed::buffer_iterator &__restrict end) {
  while ((end - it) >= 4 && is_delimiter(it[0]) && is_delimiter(it[1]) &&
         is_delimiter(it[2]) && is_delimiter(it[3])) {
    it += 4;
  }
  while (it != end && is_delimiter(*it)) {
    ++it;
  }
  if (it == end || !is_digit(*it))
    return false;

  w = 0;
  unsigned count = 0;
  while ((end - it) >= 4 && is_digit(it[0]) && is_digit(it[1]) &&
         is_digit(it[2]) && is_digit(it[3])) {
    w = w * 10 + digit(it[0]);
    w = w * 10 + digit(it[1]);
    w = w * 10 + digit(it[2]);
    w = w * 10 + digit(it[3]);
    it += 4;
    count += 4;
    if (count > MAX_DIGITS)
      throw std::runtime_error("vertex key too large");
  }

  while (it != end && is_digit(*it)) {
    if (++count > MAX_DIGITS)
      throw std::runtime_error("vertex key too large");
    w = w * 10 + digit(*it++);
  }

  return true;
}

template <typename T, typename U>
static inline void merge_sort(T &states, U &out) {
  using container_type = T::value_type;
  using iterator = container_type::iterator;
  using value_type = U::value_type;

  std::vector<std::pair<iterator, iterator>> active_ranges;

#pragma unroll 32
  for (auto &state : states) {
    if (!state.empty()) {
      active_ranges.emplace_back(state.begin(), state.end());
      logger::debug("state size", state.size());
    }
  }

  while (!active_ranges.empty()) {
    std::size_t min_idx = 0;
#pragma unroll 32
    for (std::size_t i = 1; i < active_ranges.size(); ++i) {
      if (*active_ranges[i].first < *active_ranges[min_idx].first) {
        min_idx = i;
      }
    }

    auto &[start, end] = active_ranges[min_idx];
    value_type v = *start++;
    if (out.empty() || v != out.back()) {
      out.emplace_back(std::move(v));
    }
    if (start == end) {
      if (min_idx != active_ranges.size() - 1) {
        std::swap(active_ranges[min_idx], active_ranges.back());
      }
      active_ranges.pop_back();
    }
  }
}

template <typename T> static inline void dedup(T &cont) {
  std::sort(cont.begin(), cont.end());
  cont.erase(std::unique(cont.begin(), cont.end()), cont.end());
  cont.shrink_to_fit();
}

void graph_builder::read_graph(std::pmr::vector<Edge> &edges_raw,
                               std::pmr::vector<Vertex> &vertices_raw) {
  auto stamp = std::chrono::high_resolution_clock::now();
  auto stamp_end = stamp;

  const std::size_t est_edges = 1.10 * E / T;
  const std::size_t est_vertices = 2 * est_edges;
  const std::size_t est_chunks = feed.estimated_chunks();

  std::pmr::vector<std::pmr::vector<Edge>> edges(
      T, std::pmr::vector<Edge>(memory::pool()), memory::pool());
  std::pmr::vector<std::pmr::vector<Vertex>> vertices(
      T, std::pmr::vector<Vertex>(memory::pool()), memory::pool());

#pragma unroll 32
  for (std::size_t i = 0; i < T; ++i) {
    edges[i].reserve(est_edges);
    vertices[i].reserve(est_vertices);
  }

  feed.prepare();
  for (std::size_t buff_id = 0; buff_id < est_chunks; buff_id++) {
    if (!feed) {
      logger::error("failure to read");
    }
    auto chunk = feed.read_chunk();
    pool.exec([&edges, &vertices,
               chunk = std::move(chunk)](std::uint16_t task_id) mutable {
      auto &local_edges = edges[task_id];
      auto &local_vertices = vertices[task_id];
      auto it = chunk.begin();
      auto end = chunk.end();
      while (it != end) {
        Vertex u{}, v{};
        if (!read_single_vertex(u, it, end) ||
            !read_single_vertex(v, it, end)) {
          logger::error("cannot read vertices! abort.");
        } else if (u == v) {
          logger::debug("self-loop edges are not supported");
        } else {
          local_edges.emplace_back(u, v);
          local_vertices.push_back(u);
          local_vertices.push_back(v);
        }
      }
      chunk.dispose();
    });
  }
  if (feed) {
    throw std::runtime_error("feed should be EOF by now.");
  }
  pool.join();
  feed.reclaim();

  stamp_end = std::chrono::high_resolution_clock::now();

  logger::info("Reading", logger::throughput(stamp, stamp_end, E));

  stamp = std::chrono::high_resolution_clock::now();

#pragma unroll 32
  for (std::size_t i = 0; i < T; ++i) {
    pool.exec([i, &vertices](auto) { dedup(vertices[i]); });
    pool.exec([i, &edges](auto) { dedup(edges[i]); });
  }
  pool.join();

  std::size_t edges_sz = 0;
  std::size_t vertices_sz = 0;

#pragma unroll 32
  for (std::size_t i = 0; i < T; ++i) {
    edges_sz += edges[i].size();
    vertices_sz += vertices[i].size();
  }

  edges_raw.reserve(edges_sz);
  vertices_raw.reserve(vertices_sz);

  pool.exec([&vertices_sz, &vertices, &vertices_raw](auto) {
    auto begin = std::chrono::high_resolution_clock::now();
    merge_sort(vertices, vertices_raw);
    auto end = std::chrono::high_resolution_clock::now();
    logger::info("Merging Vertices",
                 logger::throughput(begin, end, vertices_sz));
    vertices_raw.shrink_to_fit();
  });
  pool.exec([&edges_sz, &edges, &edges_raw](auto) {
    auto begin = std::chrono::high_resolution_clock::now();
    merge_sort(edges, edges_raw);
    auto end = std::chrono::high_resolution_clock::now();
    logger::info("Merging Edges", logger::throughput(begin, end, edges_sz));
    edges_raw.shrink_to_fit();
  });
  pool.join();

  stamp_end = std::chrono::high_resolution_clock::now();
  logger::info("Finished Reading Graph in",
               logger::time_diff(stamp, stamp_end));

  const std::size_t N = vertices_raw.size();
  const std::size_t M = edges_raw.size();

  if (N != V) {
    logger::warn("vertices number mismatch!", N, "should be", V);
  }
  if (M != E && M != (2 * E)) {
    logger::warn("edges number mismatch!", M, "should be", E, "or", 2 * E);
  }

  Assert(N == V, "vertices number mismatch!");
  Assert(M == E || M == 2 * E, "edges number mismatch!");

  const_cast<std::size_t &>(V) = N;
  const_cast<std::size_t &>(E) = M;
}

void graph_builder::parse_graph(
    std::pmr::vector<Vertex> &vertices_raw, std::pmr::vector<Off> &degrees,
    std::pmr::vector<std::pair<Key, Key>> &edges_key) {
  std::pmr::vector<Edge> edges_raw(memory::pool());

  read_graph(edges_raw, vertices_raw);

  auto stamp = std::chrono::high_resolution_clock::now();

  edges_key.resize(E);
  degrees.resize(V);

  std::pmr::vector<std::pmr::vector<Off>> local_degrees(
      T, std::pmr::vector<Off>(V, memory::pool()), memory::pool());

  for (std::size_t tid = 0; tid < T; ++tid) {
    std::size_t start = tid * CHUNK;
    std::size_t end = std::min(start + CHUNK, E);
    pool.exec([tid, start, end, &edges_raw, &vertices_raw, &edges_key,
               &local_degrees](std::uint16_t) {
      auto &degrees = local_degrees[tid];
      for (std::size_t i = start; i < end; ++i) {
        const auto &[u, v] = edges_raw[i];

        auto iu =
            std::lower_bound(vertices_raw.cbegin(), vertices_raw.cend(), u);
        auto u_key = Key(iu - vertices_raw.cbegin());
        edges_key[i].first = u_key;
        degrees[u_key]++;

        auto iv =
            std::lower_bound(vertices_raw.cbegin(), vertices_raw.cend(), v);
        auto v_key = Key(iv - vertices_raw.cbegin());
        edges_key[i].second = v_key;
        degrees[v_key]++;
      }
    });
  }
  pool.join();

  for (std::size_t tid = 0; tid < T; ++tid) {
    std::size_t start = tid * CHUNK;
    std::size_t end = std::min(start + CHUNK, V);
    pool.exec(
        [start, end, pool_size = T, &local_degrees, &degrees](std::uint16_t) {
          for (std::size_t t = 0; t < pool_size; ++t) {
#pragma unroll 32
            for (std::size_t i = start; i < end; ++i) {
              degrees[i] += local_degrees[t][i];
            }
          }
        });
  }
  pool.join();

  auto stamp_end = std::chrono::high_resolution_clock::now();

  logger::info("Computing Degrees", logger::throughput(stamp, stamp_end, E));
  logger::info("Finished Parsing Graph in",
               logger::time_diff(stamp, stamp_end));
}

graph graph_builder::build(bool) {
  auto begin = std::chrono::high_resolution_clock::now();
  auto stamp = begin;
  auto stamp_end = begin;

  std::pmr::vector<Vertex> vertices_raw(memory::pool());
  std::pmr::vector<Off> degrees(memory::pool());
  std::pmr::vector<std::pair<Key, Key>> edges_key(memory::pool());

  parse_graph(vertices_raw, degrees, edges_key);

  std::pmr::vector<Off> offsets(V + 1, memory::pool());
  std::vector<std::atomic<Off>> cursor(V);

#pragma unroll 32
  for (std::size_t i = 0; i < V; ++i) {
    offsets[i + 1] = offsets[i] + degrees[i];
    cursor[i].store(offsets[i], std::memory_order_relaxed);
  }

  std::pmr::vector<Key> neighs(offsets[V], memory::pool());

  for (std::size_t t = 0; t < T; ++t) {
    std::size_t start = t * CHUNK;
    std::size_t end = std::min(start + CHUNK, E);
    pool.exec([&, start, end](std::uint16_t) {
      for (std::size_t i = start; i < end; ++i) {
        auto [u, v] = edges_key[i];
        Off iu = cursor[u].fetch_add(1, std::memory_order_relaxed);
        Off iv = cursor[v].fetch_add(1, std::memory_order_relaxed);
        neighs[iu] = v;
        neighs[iv] = u;
      }
    });
  }
  pool.join();

  stamp_end = std::chrono::high_resolution_clock::now();

  logger::info("Constructing Edges", logger::throughput(stamp, stamp_end, E));

  graph G(V, E, std::move(offsets), std::move(neighs), std::move(vertices_raw));

  auto end = std::chrono::high_resolution_clock::now();
  logger::info("Graph with", V, "Vertices and", E, "Edges was Loaded in",
               logger::time_diff(begin, end));
  return G;
}

} // namespace mc
