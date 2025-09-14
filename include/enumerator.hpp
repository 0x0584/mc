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
#include "core.hpp"
#include "graph.hpp"

namespace mc {
// vertex wrapper that serves as a medium to access vertices as keys, in order
// to take advantage of std::vector (since when executing the algorithm, the
// set of vertices would remain constant) instead of how they are stored as
// std::unordered_map within the graph, it acts also as a handler of vertices
// for colouring and inducing vertex-neighbourhood
struct enumerator {

  // this is a primitive type too, same as graph::vertex so changes in the
  // implementation are required in order to avoid overhead of copying
  // instead of using references or moving the object
  using key = unsigned;
  using colour = unsigned;

  static inline const key null_key = -1u;

  class sorted_keys;

  explicit enumerator(const graph &G);

  inline std::size_t vertex_count() const { return V.size(); }

  inline graph::vertex key_to_vertex(std::size_t index) const {
    Assert(index < vertex_count());
    return V[index].first;
  }

  inline std::vector<key> neighbourhood(key v,
                                        const std::vector<key> &neighs) const {
    Assert(v < vertex_count());
    std::vector<key> new_neighs;
    new_neighs.reserve(neighs.size());
    std::copy_if(neighs.begin(), neighs.end(), std::back_inserter(new_neighs),
                 [this, v](key u) { return mtx[v].contains(u); });
    return new_neighs;
  }

  inline std::vector<key>
  neighbours(key v, std::pmr::unordered_set<key> pruned) const {
    Assert(v < vertex_count());
    std::vector<key> neighs;
    neighs.reserve(A.at(v).size());
    std::copy_if(A.at(v).begin(), A.at(v).end(), std::back_inserter(neighs),
                 [&pruned](key u) { return not pruned.contains(u); });
    return neighs;
  }

  inline std::vector<graph::vertex>
  unfold_keys(const std::vector<key> &keys) const {
    std::vector<graph::vertex> vertices(keys.size());
    std::transform(std::execution::par_unseq, keys.begin(), keys.end(),
                   vertices.begin(), [this](key v) {
                     Assert(v < vertex_count());
                     return key_to_vertex(v);
                   });
    return vertices;
  }

  void print() const;

  void draw(const std::vector<graph::vertex> &clq) const;

  sorted_keys greedy_colour_sort(std::vector<key> &&vertices) const;
  sorted_keys dsatur_colour_sort(std::vector<key> &&vertices) const;
  sorted_keys greedy_colour_sort(
      std::vector<key> &&vertices,
      lru_cache<std::vector<enumerator::key>, std::vector<enumerator::colour>>
          &cache,
      std::size_t &cache_hits) const;

  bool is_clique(const std::vector<key> &clique) const;

private:
  // FIXME: refactor this with mc::graph
  // Enumertaed vertices
  const;

  // Adjacency List for fast neighbourhood deduction
  std::vector<std::vector<key>> A;
  std::vector<std::pmr::unordered_set<key>> mtx;
};

// the keys are sorted in non-decreasing order relative to their colours
// FIXME: refactor this into a better interface
class enumerator::sorted_keys {
  std::pair<std::vector<key>, std::vector<colour>> keys_colours;

public:
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

  inline std::pair<key, colour> peek() {
    return std::make_pair(keys_colours.first.back(),
                          keys_colours.second.back());
  }

  inline std::pair<key, colour> pop() {
    auto key_colour = peek();
    keys_colours.first.pop_back();
    keys_colours.second.pop_back();
    return key_colour;
  }

  inline const std::vector<key> &keys() const { return keys_colours.first; }

  inline const std::vector<colour> &colours() const {
    return keys_colours.second;
  }

  inline std::vector<key> &keys() { return keys_colours.first; }

  inline std::vector<colour> &colours() { return keys_colours.second; }

  inline std::size_t size() const { return keys_colours.first.size(); }

  inline bool empty() const { return keys_colours.first.empty(); }

  struct iterator {

    using iterator_category = std::bidirectional_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = std::pair<key, colour>;

    inline iterator(sorted_keys *parent, std::size_t index)
        : m_parent(parent), m_index(index) {}

    [[nodiscard]] inline value_type operator*() const {
      return {m_parent->keys_colours.first[m_index],
              m_parent->keys_colours.second[m_index]};
    }

    inline iterator &operator++() {
      ++m_index;
      return *this;
    }

    [[nodiscard]] inline iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    inline iterator &operator--() {
      --m_index;
      return *this;
    }

    [[nodiscard]] inline iterator operator--(int) {
      iterator tmp = *this;
      --(*this);
      return tmp;
    }

    inline difference_type operator-(const iterator &other) const {
      return static_cast<difference_type>(m_index) -
             static_cast<difference_type>(other.m_index);
    }

    inline bool operator==(const iterator &other) const {
      return m_index == other.m_index;
    }

    inline bool operator!=(const iterator &other) const {
      return not(*this == other);
    }

  private:
    sorted_keys *m_parent;
    std::size_t m_index;

    friend class sorted_keys;
  };

  struct key_iterator {
    using iterator_category = std::bidirectional_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = key;

    inline key_iterator(iterator it) : it(it) {}

    [[nodiscard]] inline value_type operator*() const { return (*it).first; }

    inline key_iterator &operator++() {
      ++it;
      return *this;
    }

    [[nodiscard]] inline key_iterator operator++(int) {
      key_iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    inline key_iterator &operator--() {
      --it;
      return *this;
    }

    [[nodiscard]] inline key_iterator operator--(int) {
      key_iterator tmp = *this;
      --(*this);
      return tmp;
    }

    inline difference_type operator-(const key_iterator &other) const {
      return it - other.it;
    }

    inline bool operator==(const key_iterator &other) const {
      return it == other.it;
    }

    inline bool operator!=(const key_iterator &other) const {
      return not(*this == other);
    }

  private:
    iterator it;

    friend class sorted_keys;
  };

  iterator begin() { return iterator(this, 0); }
  iterator end() { return iterator(this, size()); }
};
} // namespace mc

#endif // ENUMERATOR_HPP
