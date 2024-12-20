// graph.hpp
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

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <cassert>

#include <ostream>
#include <future>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory_resource>

#include "log.hpp"
#include "input.hpp"

namespace mc {
struct graph {
  // this is a primitive type, ins case of a change in the implementation
  // since some parts of the code should be updated to avoid overhead of
  // copying the  objects rather than either referencing them or moving them
  using vertex = unsigned;

  using neighbours_set = std::pmr::unordered_set<vertex>;
  using adjacency_map = std::pmr::unordered_map<vertex, neighbours_set>;

  friend struct enumerator;
  friend struct graph_builder;

  graph(const graph &) = delete;
  graph(graph &&) = default;
  explicit graph(std::pmr::monotonic_buffer_resource &vertices_pool)
      : A(&vertices_pool) {}
  graph &operator=(graph &) = delete;
  graph &operator=(graph &&G) = default;

  inline const neighbours_set &neighbours(vertex v) const { return A.at(v); }
  inline const adjacency_map &adjacency() const { return A; }
  inline bool directed() { return not undirected; }

  void print() const;

private:
  bool undirected = false;
  std::size_t edge_count = 0;
  adjacency_map A;
};

struct graph_builder {
  graph_builder(graph_builder &&) = delete;
  graph_builder(const graph_builder &) = delete;

  explicit graph_builder(input &in,
                         std::pmr::monotonic_buffer_resource &vertices_pool,
                         std::pmr::monotonic_buffer_resource &edges_pool,
                         bool undirected = args::undirected)
      : feed(in), G(vertices_pool), edges_pool(edges_pool) {
    G.undirected = undirected;
    Q.reserve(feed.estimate_chunks());
  }

  graph_builder &operator=(const graph_builder &) = delete;
  graph_builder &operator=(graph_builder &&) = delete;

  bool read_single_vertex(graph::vertex &w, std::string::iterator &it,
                          std::string::iterator end);

  graph build();

private:
  std::mutex graph_mtx;
  std::vector<std::future<void>> Q;
  mc::feed feed;
  graph G;
  std::pmr::monotonic_buffer_resource &edges_pool;
};
} // namespace mc
#endif // GRAPH_HPP