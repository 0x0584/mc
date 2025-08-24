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
#include "pq.hpp"
#include <algorithm>
#include <memory_resource>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace mc {
struct vertex {
  using value_type = std::uint64_t;

  explicit vertex(value_type val) : val(val) {}

  vertex() = default;
  vertex(const vertex &) = default;
  vertex(vertex &&) = default;

  operator value_type() const { return val; }
  vertex &operator=(const vertex &) = default;
  vertex &operator=(vertex &&) = default;
  bool operator==(const vertex &) const = default;

private:
  value_type val = -1u;
};

struct colour {
  using value_type = std::uint64_t;

  explicit colour(value_type val) : val(val) {}

  colour() = default;
  colour(const colour &) = default;
  colour(colour &&) = default;

  bool operator==(const colour &) const = default;
  colour &operator++() {
    val += 1;
    return *this;
  }

  operator value_type() const { return val; }
  colour &operator=(const colour &) = default;
  colour &operator=(colour &&) = default;

  bool is_nil() const { return val == -1ul; }

private:
  value_type val = -1ul;
};
} // namespace mc

namespace std {
template <> struct hash<mc::vertex> {
  [[nodiscard]] inline std::size_t operator()(const mc::vertex &v) const {
    return hash_func(v);
  }

private:
  std::hash<mc::vertex::value_type> hash_func;
};

template <> struct hash<mc::colour> {
  [[nodiscard]] inline std::size_t operator()(const mc::colour &v) const {
    return hash_func(v);
  }

private:
  std::hash<mc::colour::value_type> hash_func;
};

} // namespace std

namespace mc {
// FIXME: refactor the graph class
struct graph {
  friend struct graph_builder;

  // this is a primitive type, ins case of a change in the implementation
  // since some parts of the code should be updated to avoid overhead of
  // copying the  objects rather than either referencing them or moving them

  using key = std::uint64_t;

  graph() = default;
  graph(const graph &) = delete;
  graph(graph &&) = default;

  ~graph() { logger::debug("~graph()"); }

  graph &operator=(graph &) = delete;
  graph &operator=(graph &&) = default;

  inline std::size_t vertex_count() const { return V.size(); }
  inline std::size_t edge_count() const { return _edge_count; }

  inline vertex to_vertex(key k) const { return V[k]; }

  inline std::pmr::vector<vertex>
  to_vertex(const std::pmr::vector<key> &keys) const {
    std::pmr::vector<vertex> verts(memory::pool());
    verts.reserve(keys.size());
    for (key k : keys) {
      verts.emplace_back(V[k]);
    }
    return verts;
  }

  inline key to_key(const vertex &v) const { return vertex_to_key.at(v); }

  inline const std::pmr::unordered_set<key> &neighbours(key k) const {
    return A[k];
  }

  inline std::pmr::vector<key>
  neighbours(key v, const std::pmr::vector<key> &neighs) const {
    std::pmr::vector<key> new_neighs(memory::pool());
    new_neighs.reserve(neighs.size());
    for (key u : neighs) {
      if (A[v].contains(u)) {
        new_neighs.emplace_back(u);
      }
    }
    return new_neighs;
  }

  bool is_clique(const std::pmr::vector<key> &clique) const {
    for (key u : clique) {
      for (key v : clique) {
        if (u == v) {
          continue;
        } else if (!A[u].contains(v)) {
          return false;
        }
      }
    }
    return true;
  }

  std::pair<std::pmr::vector<key>, std::pmr::vector<colour>>
  colour_sort(const std::pmr::vector<key> &neighs) const {
    make_scope_timer(colour_sort_timer);

    // TODO: use thread local vector instead
    thread_local std::pmr::vector<colour> colours(vertex_count(),
                                                  memory::pool());
    thread_local std::pmr::vector<std::size_t> saturation(vertex_count(),
                                                          memory::pool());
    thread_local auto sort_keys = [this](key u, key v) {
      return saturation[u] != saturation[v] ? saturation[u] < saturation[v]
                                            : A[u].size() < A[v].size();
    };
    thread_local utils::pq<std::size_t, decltype(sort_keys)> q(sort_keys);

    {
      make_scope_timer(setup_timer);
      std::fill(colours.begin(), colours.end(), colour());
      std::fill(saturation.begin(), saturation.end(), 0);
    }
    {
      make_scope_timer(queue_emplace_timer);
      for (key u : neighs) {
        q.emplace(u);
      }
    }

    std::pmr::vector<std::pair<key, colour>> sorted(neighs.size(),
                                                    memory::pool());

    while (!q.empty()) {
      key u = q.top();

      std::pmr::vector<bool> used_colours(A[u].size() + 1, memory::pool());
      for (key v : A[u]) {
        if (!colours[v].is_nil()) {
          used_colours[colours[v]] = true;
        }
      }

      colour col(0);
      while (col < A[u].size() && used_colours[col]) {
        ++col;
      }
      colours[u] = col;

      sorted.emplace_back(u, col);

      // for (key v : A[u]) {
      //   if (colours[v].is_nil()) {
      //     std::size_t count = 0;
      //     for (key w : A[v]) {
      //       if (!colours[w].is_nil()) {
      //         count++;
      //       }
      //     }
      //     saturation[v] = count;
      //   }
      // }

      for (key v : A[u]) {
        if (colours[v].is_nil()) {
          bool is_new_color = true;
          for (key w : A[v]) {
            if (!colours[w].is_nil() && colours[w] == colours[u]) {
              is_new_color = false;
              break;
            }
          }
          if (is_new_color) {
            saturation[v]++;
          }
        }
      }

      {
        make_scope_timer(pop_timer);
        q.pop();
      }
    }

    {
      make_scope_timer(sort_timer);
      std::sort(sorted.begin(), sorted.end(), [](const auto &a, const auto &b) {
        return a.second < b.second;
      });
    }
    std::pmr::vector<key> keys(memory::pool());
    keys.reserve(sorted.size());
    std::pmr::vector<colour> cols(memory::pool());
    cols.reserve(sorted.size());

    {
      make_scope_timer(copy_timer);
      for (auto &[key, colour] : sorted) {
        keys.emplace_back(std::move(key));
        cols.emplace_back(std::move(colour));
      }
    }

    return std::make_pair(std::move(keys), std::move(cols));
  }

  void print() const;

  /*
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
   */
private:
  std::size_t _edge_count = 0;

  std::pmr::unordered_map<vertex, key> vertex_to_key;
  std::pmr::vector<vertex> V;
  std::pmr::vector<std::pmr::unordered_set<key>> A;
};

struct graph_builder {
  graph_builder(graph_builder &&) = delete;
  graph_builder(const graph_builder &) = delete;

  explicit graph_builder(input &in) : feed(in) {}

  ~graph_builder() { logger::debug("~graph_builder()"); }

  graph_builder &operator=(const graph_builder &) = delete;
  graph_builder &operator=(graph_builder &&) = delete;

  bool read_single_vertex(vertex &w, std::string::iterator &it,
                          std::string::iterator end);

  graph build(bool undirected);

private:
  mc::feed feed;
};
} // namespace mc
#endif // GRAPH_HPP
