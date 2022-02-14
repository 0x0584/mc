#pragma once

#include <assert.h>

#include <iostream>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <utility>
#include <functional>
#include <vector>
#include <algorithm>

template<typename T> inline void hash_combine(std::size_t &seed, T const &key) {
  std::hash<T> hasher;
  seed ^= hasher(key) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
  template<typename T1, typename T2> struct hash<std::pair<T1, T2>> {
    std::size_t operator()(std::pair<T1, T2> const &p) const {
      std::size_t seed1(0);
      ::hash_combine(seed1, p.first);
      ::hash_combine(seed1, p.second);

      std::size_t seed2(0);
      ::hash_combine(seed2, p.second);
      ::hash_combine(seed2, p.first);

      return std::min(seed1, seed2);
    }
  };
}

// template<typename Iterator, class Function>
// void parallel_for(const Iterator& first, const Iterator& last, Function &&f,
//                   const int nthreads = std::thread::hardware_concurrency(), const int threshold = 1) {
//   const unsigned int group = std::max(std::max(1, std::abs(threshold)), (last-first)/std::abs(nthreads));
//   std::vector<std::thread> threads;
//   threads.reserve(nthreads);
//   Iterator it = first;
//   for (; last - it > group; it += group) {
//     threads.push_back(std::thread([=, &f] {std::for_each(it, std::min(it + group, last), f);}));
//   }
//   threads.push_back(std::thread([=, &f] {std::for_each(it, last, f);}));
//   std::for_each(threads.begin(), threads.end(),
//                 [](std::thread& x) {x.join();});
// }

using vertex_t = int;

struct edge {
  edge(const edge &e) = default;
  edge(vertex_t u, vertex_t v) : s(u), t(v) {}

  edge residual() const { return {t, s}; }
  vertex_t from() const { return s; }
  vertex_t to() const { return t; }

  edge &operator=(const edge &) = delete;

private:
  vertex_t s, t;
};

struct graph {
  using adjacency_matrix_t = std::vector<std::vector<bool>>;
  using adjacency_list_t = std::unordered_map<vertex_t, std::vector<vertex_t>>;
  using color_t = int;

  graph() = default;
  graph(const graph &) = delete;
  graph(const std::vector<vertex_t> &V,
        const std::vector<edge> &E,
        bool directed = true);

  std::vector<vertex_t> vertices() const;
  const std::vector<vertex_t> &neighbors(vertex_t u) const;
  std::size_t num_vertices() const { return edges.size(); }
  std::size_t degree(vertex_t u) const;
  adjacency_matrix_t adjacency_matrix() const;
  const adjacency_list_t &adjacency_list() const;

  bool add_vertex(vertex_t v);
  bool add_edge(const edge &e);
  void add_edge_undirected(const edge &e);

  bool remove_vertex(vertex_t v);
  bool remove_edge(const edge &e);
  void remove_edge_undirected(const edge &e);

  void print() const;

  graph &operator=(const graph &) = delete;

  static graph generate();

  std::vector<vertex_t> fahle_max_clique() const;
  std::vector<vertex_t> tomita_max_clique() const;

  std::pair<std::vector<vertex_t>, std::vector<color_t>> color_welsh_powell() const;

private:
  std::unordered_map<vertex_t, std::vector<vertex_t>> edges;
};
