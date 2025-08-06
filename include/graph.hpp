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

#include <memory_resource>
#include <unordered_map>
#include <unordered_set>

#include "input.hpp"

namespace mc {
struct graph {
  // this is a primitive type, ins case of a change in the implementation
  // since some parts of the code should be updated to avoid overhead of
  // copying the  objects rather than either referencing them or moving them
  using vertex = unsigned;

  static inline const vertex nil_vertex = -1u;

  using neighbours_set = std::pmr::unordered_set<vertex>;
  using adjacency_map = std::pmr::unordered_map<vertex, neighbours_set>;

  friend struct enumerator;
  friend struct graph_builder;

  explicit inline graph(std::pmr::monotonic_buffer_resource &vertices_pool,
                        std::pmr::monotonic_buffer_resource &edges_pool)
      : vertices_pool(vertices_pool), edges_pool(edges_pool),
        A(&vertices_pool) {}

  ~graph() { logger::debug("~graph()"); }

  graph(const graph &) = delete;
  graph(graph &&) = default;

  graph &operator=(graph &) = delete;
  graph &operator=(graph &&) = delete;

  inline const neighbours_set &neighbours(vertex v) const { return A.at(v); }
  inline const adjacency_map &adjacency() const { return A; }
  inline bool directed() { return not undirected; }
  inline std::size_t vertex_count() const { return A.size(); }
  inline std::size_t edge_count() const { return _edge_count; }
  inline const adjacency_map &adjacency_list() const { return A; }
  bool add_edge_undirected(vertex u, vertex v, std::size_t edge_set_size = 0);

  void print() const;

private:
  bool undirected = false;
  std::size_t _edge_count = 0;

  std::pmr::monotonic_buffer_resource &vertices_pool;
  std::pmr::monotonic_buffer_resource &edges_pool;

  adjacency_map A;
};

struct graph_builder {
  graph_builder(graph_builder &&) = delete;
  graph_builder(const graph_builder &) = delete;

  explicit graph_builder(input &in) : feed(in) {}

  ~graph_builder() { logger::debug("~graph_builder()"); }

  graph_builder &operator=(const graph_builder &) = delete;
  graph_builder &operator=(graph_builder &&) = delete;

  bool read_single_vertex(graph::vertex &w, std::string::iterator &it,
                          std::string::iterator end);

  graph build(std::pmr::monotonic_buffer_resource &vertices_pool,
              std::pmr::monotonic_buffer_resource &edges_pool);

private:
  mc::feed feed;
};
} // namespace mc
#endif // GRAPH_HPP
