// enumerator.hpp
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

#ifndef ENUMERATOR_HPP
#define ENUMERATOR_HPP

#include <algorithm>
#include <execution>
#include <shared_mutex>

#include "cache.hpp"
#include "graph.hpp"

namespace mc {
// vertex wrapper that serves as a medium to access vertices as keys, in order
// to take advantage of std::vector (since when executing the algorithm, the
// set of vertices would remain constant) instead of how they are stored as
// std::unordered_map within the graph, it acts also as a handler of vertices
// for colouring and inducing vertex-neighbourhood
struct enumerator {
  using adjacency_vector =
      std::vector<std::pair<graph::vertex, graph::neighbours_set>>;

  // this is a primitive type too, same as graph::vertex so changes in the
  // implementation are required in order to avoid overhead of copying
  // instead of using references or moving the object
  using key = unsigned;
  using colour = unsigned;

  // the keys are sorted in non-decreasing order relative to their colours
  struct sorted_keys { // FIXME: refactor this into a better interface
    friend enumerator;

    sorted_keys() = default;

    explicit inline sorted_keys(
        std::pair<std::vector<key>, std::vector<colour>> sorted)
        : keys_colours(std::move(sorted)) {}

    inline colour highest_colour() const { return keys_colours.second.back(); }

    inline key key_with_highest_colour() const {
      return keys_colours.first.back();
    }

    inline key pop_key_with_highest_colour() {
      key k = keys_colours.first.back();
      keys_colours.first.pop_back();
      keys_colours.second.pop_back();
      return k;
    }

    inline const std::vector<key> &keys() const { return keys_colours.first; }

    inline const std::vector<colour> &colours() const {
      return keys_colours.second;
    }

    inline std::size_t size() const { return keys_colours.first.size(); }

    inline bool empty() const { return keys_colours.first.empty(); }

  private:
    std::pair<std::vector<key>, std::vector<colour>> keys_colours;
  };

  explicit enumerator(graph &G);

  inline std::size_t vertex_count() const { return V.size(); }

  inline std::size_t get_cache_hits() { return cache_hits; }

  inline void reset_cache_hits() { cache_hits = 0; }

  inline graph::vertex key_to_vertex(std::size_t index) const {
    assert(index < vertex_count());
    return V[index].first;
  }

  inline std::vector<key> neighbourhood(key v,
                                        const std::vector<key> &neighs) const {
    assert(v < vertex_count());
    std::vector<key> new_neighs;
    new_neighs.reserve(neighs.size());
    std::copy_if(neighs.begin(), neighs.end(), std::back_inserter(new_neighs),
                 // neighbours of both vertices u and v
                 [this, v](key u) { return B[v][u]; });
    return new_neighs;
  }

  inline std::vector<key>
  neighbours(key v, std::pmr::unordered_set<key> &pruned) const {
    assert(v < vertex_count());
    std::vector<key> neighs;
    neighs.reserve(A.at(v).size());
    std::copy_if(A.at(v).begin(), A.at(v).end(), std::back_inserter(neighs),
                 [&pruned](key u) { return not pruned.count(u); });
    return neighs;
  }

  inline std::vector<graph::vertex>
  unfold_keys(const std::vector<key> &keys) const {
    std::vector<graph::vertex> vertices(keys.size());
    std::transform(std::execution::par_unseq, keys.begin(), keys.end(),
                   vertices.begin(), [this](key v) {
                     assert(v < vertex_count());
                     return key_to_vertex(v);
                   });
    return vertices;
  }

  void print() const;

  void draw(const std::vector<graph::vertex> &clq) const;

  sorted_keys greedy_colour_sort(std::vector<key> &&vertices) const;

  bool is_clique(const std::vector<key> &clique) const;

private:
  void cache_hit_progress() const;

  mutable std::shared_mutex cache_mtx;
  mutable lru_cache<std::vector<key>, std::vector<colour>> cache;
  mutable std::atomic_size_t cache_hits{0};

  adjacency_vector V; // Enumertaed vertices

  // Adjacency List for fast neighbourhood deduction
  std::vector<std::vector<key>> A;
  // Adjacency Matrix for fast edge probing
  std::vector<std::vector<bool>> B;
};
} // namespace mc

#endif // ENUMERATOR_HPP
