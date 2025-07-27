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

#include "graph.hpp"
#include <algorithm>
#include <execution>
#include <functional>
#include <list>
#include <set>
#include <shared_mutex>

namespace mc {
// vertex wrapper that serves as a medium to access vertices as keys, in order
// to take advantage of std::vector (since when executing the algorithm, the
// set of vertices would remain constant) instead of how they are stored as
// std::unordered_map within the graph, it acts also as a handler of vertices
// for colouring and inducing vertex-neighbourhood
struct enumerator {
  using adjacency_vector =
      std::vector<std::pair<graph::vertex, graph::neighbours_set>>;
  // this is a premitive type too, same as graph::vertex so changes in the
  // implementation are required in order to avoid overhead of copying
  // instead of using references or moving the object
  using key = unsigned;
  using colour = unsigned;

  // the keys are sorted in non-decreasing order relative to their colours
  struct sorted_keys {
    friend enumerator;

    sorted_keys() {};
    explicit sorted_keys(
        std::pair<std::vector<key>, std::vector<colour>> sorted)
        : keys_colours(std::move(sorted)) {}

    colour highest_colour() const { return keys_colours.second.back(); }

    key key_with_highest_colour() const { return keys_colours.first.back(); }
    key pop_key_with_highest_colour() {
      key k = keys_colours.first.back();
      keys_colours.first.pop_back();
      keys_colours.second.pop_back();
      return k;
    }

    const std::vector<key> &keys() const { return keys_colours.first; }

    std::size_t size() const { return keys_colours.first.size(); }

    bool empty() const { return keys_colours.first.empty(); }

    std::string to_string() const {
      std::ostringstream oss;
      oss << "{ ";
      for (auto i = 0ul; i < size(); ++i) {
        oss << "'" << keys_colours.first[i] << "':" << keys_colours.second[i]
            << " ";
      }
      oss << "}";
      return oss.str();
    }

  private:
    std::pair<std::vector<key>, std::vector<colour>> keys_colours;
  };

  explicit enumerator(graph G);

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

  void cache_hit_progress() const;

  sorted_keys greedy_colour_sort(std::vector<key> &&vertices) const;

  bool is_clique(const std::vector<key> &clique) const;

  inline std::size_t vertex_count() const { return V.size(); }

  std::size_t get_cache_hits() { return cache_hits; }
  void reset_cache_hits() { cache_hits = 0; }

  // private:
  adjacency_vector V; // Enumertaed vertices

  template <typename T> static void print_vector(const std::vector<T> &v) {
    for (T i : v) {
      std::cout << i << " ";
    }
    std::cout << "\n";
  }

  struct colouring_cache {
    template <typename T> struct vector_hash {
      std::size_t
      operator()(const std::reference_wrapper<std::vector<T>> &v_ref) const {
        std::size_t seed = v_ref.get().size();
        std::vector<T> tmp = v_ref.get();
        std::sort(tmp.begin(), tmp.end());
        for (T e : tmp) {
          seed ^= hash_func(e) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
      }

    private:
      std::hash<T> hash_func;
    };

    template <typename T> struct vector_equal {
      bool operator()(const std::reference_wrapper<std::vector<T>> &a,
                      const std::reference_wrapper<std::vector<T>> &b) const {
        std::vector<T> atmp = a.get();
        std::sort(atmp.begin(), atmp.end());
        std::vector<T> btmp = b.get();
        std::sort(btmp.begin(), btmp.end());
        return atmp == btmp;
      }
    };

    explicit colouring_cache(std::size_t capacity)
        : _capacity(capacity), keys(&buff), values(&buff) {
      assert(capacity > 0);
      keys.reserve(capacity);
    }

    std::optional<std::pair<std::vector<key>, std::vector<colour>>>
    get(std::vector<key> &k) {
      std::unique_lock<std::mutex> lock(mtx);
      std::optional<std::pair<std::vector<key>, std::vector<colour>>> value;
      if (keys.contains(k)) {
        std::pmr::list<
            std::pair<std::vector<key>, std::vector<colour>>>::iterator it =
            keys[k];
        cache_hit(it);
        value.emplace(*it);
      }
      return value;
    }

    std::pair<std::vector<key>, std::vector<colour>>
    set(std::vector<key> &k, std::vector<colour> &v) {
      std::unique_lock<std::mutex> lock(mtx);
      std::pmr::list<std::pair<std::vector<key>, std::vector<colour>>>::iterator
          it;
      if (keys.contains(k)) {
        it = keys[k];
        it->second = std::move(v);
        cache_hit(it);
      } else {
        if (values.size() == _capacity) {
          it = values.end();
          --it;
          keys.erase(it->first);
          it->first = std::move(k);
          it->second = std::move(v);
          cache_hit(it);
        } else {
          it = values.emplace(values.begin(), std::move(k), std::move(v));
          if (values.size() == _capacity) {
            log::info(COL_MAGENTA, "cache is full!");
          }
        }
        keys.emplace(it->first, it);
      }
      return *it;
    }

    std::size_t size() const {
      std::unique_lock<std::mutex> lock(mtx);
      return values.size();
    }

    std::size_t capacity() const { return _capacity; }

    bool full() const { return size() == _capacity; }

  private:
    inline void
    cache_hit(std::pmr::list<
              std::pair<std::vector<key>, std::vector<colour>>>::iterator it) {
      values.splice(values.begin(), values, it);
    }

    mutable std::mutex mtx;
    const std::size_t _capacity;

    std::pmr::monotonic_buffer_resource buff;
    std::pmr::unordered_map<
        std::reference_wrapper<std::vector<key>>,
        std::pmr::list<
            std::pair<std::vector<key>, std::vector<colour>>>::iterator,
        vector_hash<key>, vector_equal<key>>
        keys;
    std::pmr::list<std::pair<std::vector<key>, std::vector<colour>>> values;
  };

  mutable std::shared_mutex cache_mtx;
  mutable colouring_cache cache;

  mutable std::atomic_size_t cache_hits{0};

  template <typename T> using vector_2d = std::vector<std::vector<T>>;
  vector_2d<key> A;  // Adjacency List for fast neighbourhood deduction
  vector_2d<bool> B; // Adjacency Matrix for fast edge probing
};
} // namespace mc

#endif // ENUMERATOR_HPP
