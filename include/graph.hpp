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
struct graph {
  friend struct graph_builder;

  typedef uint32_t vertex;

  graph() = delete;
  graph(const graph &) = delete;
  graph(graph &&) = default;

  ~graph() { logger::debug("~graph()"); }

  graph &operator=(graph &) = delete;
  graph &operator=(graph &&) = delete;

  inline std::size_t vertex_count() const { return n_vertices; }

  inline std::size_t edge_count() const { return n_edges; }

  // FIXME: benchmark storing them separately
  inline degree degree(const key u) const {
    return static_cast<mc::degree>(offsets[u + 1] - offsets[u]);
  }

  inline const vertex &to_vertex(key u) const { return key_to_vertex[u]; }

  inline std::pmr::vector<vertex>
  to_vertex(const std::pmr::vector<key> &keys) const {
    std::pmr::vector<vertex> vertices(memory::pool());
    vertices.reserve(keys.size());
    for (const auto &u : keys) {
      vertices.emplace_back(key_to_vertex[u]);
    }
    return vertices;
  }

  inline std::pair<const key *, const key *> neighbours(key u) const {
    const offset first = offsets[u];
    const offset last = offsets[u + 1];
    return {neighs.data() + first, neighs.data() + last};
  }

  inline std::pmr::vector<key>
  neighbourhood(key u, const std::pmr::vector<key> &keys) const {
    std::pmr::vector<key> neighs(memory::pool());
    neighs.reserve(degree(u));
    auto [first, last] = neighbours(u);
    for (const auto &v : keys) {
      if (std::binary_search(first, last, v)) {
        neighs.emplace_back(v);
      }
    }
    return neighs;
  }

  inline bool is_clique(const std::pmr::vector<key> &keys) const {
    const size_t k = keys.size();
    if (k <= 1)
      return true;

    std::atomic<bool> failed{false};
    thread::pool pool(args::num_threads);

    for (size_t i = 0; i < k; ++i) {
      pool.submit([&, i](auto) {
        if (failed.load(std::memory_order_acquire))
          return;

        const key u = keys[i];

        if (degree(u) < k - 1) {
          failed.store(true, std::memory_order_release);
          return;
        }

        auto [first, last] = neighbours(u);

        for (size_t j = i + 1; j < k; ++j) {
          if (failed.load(std::memory_order_acquire))
            return;

          const key v = keys[j];
          if (!std::binary_search(first, last, v)) {
            failed.store(true, std::memory_order_release);
            return;
          }
        }
      });
    }

    pool.join();
    return !failed.load(std::memory_order_acquire);
  }

  void print() const;

private:
  graph(std::size_t n_vertices, std::size_t n_edges,
        std::pmr::vector<offset> offsets, std::pmr::vector<key> neighs,
        std::pmr::vector<vertex> key_to_vertex)
      : n_vertices(n_vertices), n_edges(n_edges), offsets(std::move(offsets)),
        neighs(std::move(neighs)), key_to_vertex(std::move(key_to_vertex)) {
    verify_invariants();
  }

  const std::size_t n_vertices;
  const std::size_t n_edges;
  const std::pmr::vector<offset> offsets;
  const std::pmr::vector<key> neighs;
  const std::pmr::vector<vertex> key_to_vertex;

  void verify_invariants() {
#ifndef NDEBUG
    const std::size_t V = vertex_count();
    const std::size_t E = edge_count();

    // The internal arrays are private, so we'll access them via public methods,
    // or assume friendly access for this verification utility.
    // Assuming internal access:
    // const auto& offsets = offsets;
    // const auto& neighs = neighs;

    logger::info("Graph Verification: Starting comprehensive check...");

    // --- 1. Size Invariants Check ---
    // The total size of the adjacency list must equal the stated edge count.
    Assert(offsets[V] == 2 * E,
           "VERIFY: Total length of adjacency lists (offsets[V]=" +
               std::to_string(offsets[V]) + ") != stated edge count (E=" +
               std::to_string(2 * E) + "). CSC structure is broken.");

    // --- 2. Iterative Adjacency List Check (Duplicates & Sorting) ---
    std::size_t verified_edge_count = 0;

    for (mc::key u = 0; u < V; ++u) {
      // Retrieve the neighbors range using the public method
      auto [first, last] = neighbours(u);
      const auto list_size = std::distance(first, last);

      // a. Degree/Offset Consistency
      Assert(degree(u) == list_size,
             "VERIFY: Degree method (" + std::to_string(degree(u)) +
                 ") != Adjacency list size (" + std::to_string(list_size) +
                 ") for vertex " + std::to_string(u));

      // b. Self-Loop Check (Simple Graph Invariant)
      if (list_size > 0) {
        Assert(*first != u, "VERIFY: Self-loop found for vertex " +
                                std::to_string(u) + " at start of list.");
      }

      // c. Sorted & Uniqueness Check (CRITICAL for the d=61 error)
      mc::key prev_v = 0;
      bool is_first = true;

      for (auto it = first; it != last; ++it) {
        const mc::key v = *it;

        // i. Key Validity Check
        Assert(v < V, "VERIFY: Invalid neighbor key v=" + std::to_string(v) +
                          " found in u=" + std::to_string(u) + "'s list.");

        // ii. Sorted Check (Check the sort function worked)
        if (!is_first) {
          Assert(prev_v <= v, "VERIFY: Adjacency list is NOT sorted for u=" +
                                  std::to_string(u) +
                                  ". prev=" + std::to_string(prev_v) +
                                  ", curr=" + std::to_string(v));
        }

        // iii. Uniqueness Check (Detects the d=61/multi-edge source)
        if (!is_first) {
          // If this assertion fails, the graph builder failed to deduplicate.
          Assert(
              prev_v != v,
              "VERIFY: Duplicate neighbor (multi-edge) v=" + std::to_string(v) +
                  " found in u=" + std::to_string(u) +
                  "'s list. This is the likely cause of the d > R_size error.");
        }

        prev_v = v;
        is_first = false;
      }

      verified_edge_count += std::size_t(list_size);
    }

    // Final check for total edge count consistency.
    Assert(verified_edge_count == 2 * E,
           "VERIFY: Total verified edges (" +
               std::to_string(verified_edge_count) +
               ") does not match stated edge count (E=" + std::to_string(E) +
               ").");

    logger::info("Graph Verification: All checks passed. Graph structure is "
                 "confirmed valid.");
#endif
  }
};

struct graph_builder {
  using Vertex = graph::vertex;
  using Key = key;
  using Off = offset;
  using Edge = std::pair<Vertex, Vertex>;

  graph_builder(graph_builder &&) = delete;
  graph_builder(const graph_builder &) = delete;

  explicit graph_builder(input &in)
      : feed(in), T(args::num_threads), V(feed.num_vertices()),
        E(feed.num_edges()), pool(T) {}

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

  thread::pool pool;
};
} // namespace mc
#endif // GRAPH_HPP
