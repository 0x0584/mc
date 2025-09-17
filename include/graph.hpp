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

#include "input.hpp"
#include <algorithm>

namespace mc {
// FIXME: refactor the graph class
// graphs are stored as Column Sparse Representation for optimal storage and
// manipulation.  original vertices are sorted and stored in `key_to_vertex'
// which serves as a mapping between keys and vertices.  a key is the internal
// manner in which the graph manipulates the vertices and edges.  the edges are
// stored in `neighs', all graphs are considered undirected.  the neighbours of
// a particular vertex v with key k are the subarray indicated by `offsets' as
// folows `N(v)=[offsets[k]..offsets[k+1]]'.
struct graph {
  friend struct graph_builder;

  using vertex = std::uint32_t;
  using key = std::uint32_t;
  using offset = std::uint32_t;
  using colour = std::uint8_t;

  struct vertex_subset_index {
    explicit vertex_subset_index(const std::pmr::vector<key> &vs)
        : vs_sorted(vs) {
      std::sort(vs_sorted.begin(), vs_sorted.end());
      // vs_sorted.erase(std::unique(vs_sorted.begin(), vs_sorted.end()),
      //                 vs_sorted.end());
    }

    bool contains(key u) const {
      return std::binary_search(vs_sorted.cbegin(), vs_sorted.cend(), u);
    }

    std::size_t index(key u) const {
      auto it = std::lower_bound(vs_sorted.cbegin(), vs_sorted.cend(), u);
      return static_cast<std::size_t>(it - vs_sorted.cbegin());
    }

    const std::pmr::vector<key> &vertices() const { return vs_sorted; }

    std::size_t size() const { return vs_sorted.size(); }

  private:
    std::pmr::vector<key> vs_sorted;
  };

  graph() = delete;
  graph(const graph &) = delete;
  graph(graph &&) = default;

  ~graph() { logger::debug("~graph()"); }

  graph &operator=(graph &) = delete;
  graph &operator=(graph &&) = default;

  inline std::size_t vertex_count() const { return n_vertices; }

  inline std::size_t edge_count() const { return n_edges; }

  inline offset degree(const key u) const {
    return offsets[u + 1] - offsets[u];
  }

  inline const vertex &to_vertex(key u) const { return key_to_vertex[u]; }

  inline std::pair<std::pmr::vector<key>::const_iterator,
                   std::pmr::vector<key>::const_iterator>
  neighbours(key u) const {
    const offset first = offsets[u];
    const offset last = offsets[u + 1];
    return {neighs.cbegin() + first, neighs.cbegin() + last};
  }

  bool is_clique(const std::pmr::vector<key> &vs) const {
    vertex_subset_index idx(vs);
    auto const &S = idx.vertices();

    std::size_t k = idx.size();
    if (k < 2)
      return true;

    for (key u : S) {
      if (degree(u) < k - 1)
        return false;
    }

    for (key u : S) {
      std::size_t cnt = 0;
      auto [b, e] = neighbours(u);
      for (auto it = b; it != e; ++it) {
        if (idx.contains(*it))
          ++cnt;
      }
      if (cnt != k - 1)
        return false;
    }
    return true;
  }

  void print() const;

private:
  graph(std::size_t n_vertices, std::size_t n_edges,
        std::pmr::vector<offset> offsets, std::pmr::vector<key> neighs,
        std::pmr::vector<vertex> key_to_vertex)
      : n_vertices(n_vertices), n_edges(n_edges), offsets(std::move(offsets)),
        neighs(std::move(neighs)), key_to_vertex(std::move(key_to_vertex)) {}

  std::size_t n_vertices;
  std::size_t n_edges;

  std::pmr::vector<offset> offsets;
  std::pmr::vector<key> neighs;
  std::pmr::vector<vertex> key_to_vertex;
};

struct graph_builder {
  using Vertex = graph::vertex;
  using Key = graph::key;
  using Off = graph::offset;
  using Edge = std::pair<Vertex, Vertex>;

  graph_builder(graph_builder &&) = delete;
  graph_builder(const graph_builder &) = delete;

  explicit inline graph_builder(input_source &in)
      : feed(in), T(args::num_threads), V(feed.num_vertices()),
        E(feed.num_edges()), CHUNK((E + T - 1) / T), pool(T) {}

  ~graph_builder() { logger::debug("~graph_builder()"); }

  graph_builder &operator=(const graph_builder &) = delete;
  graph_builder &operator=(graph_builder &&) = delete;

  graph build(bool undirected);

private:
  void read_graph(std::pmr::vector<Edge> &edges_raw,
                  std::pmr::vector<Vertex> &vertices_raw);
  void parse_graph(std::pmr::vector<Vertex> &vertices_raw,
                   std::pmr::vector<Off> &degrees,
                   std::pmr::vector<std::pair<Key, Key>> &edges_key);

  mc::feed feed;

  const std::size_t T;
  const std::size_t V;
  const std::size_t E;
  const std::size_t CHUNK;

  thread::pool pool;
};
} // namespace mc
#endif // GRAPH_HPP
