#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "utils.hpp"

namespace mc {
struct graph {
  friend struct reader;

  using vertex = size_type;
  using vertex_set = std::vector<vertex>;

  graph() = default;
  graph(const graph &) = delete;
  graph(graph &&) = default;

  size_type num_vertices() const { return num_v; }
  size_type num_edges() const { return num_e; }

  bool adjacent(const vertex v, const vertex u) const {
    assert(v < num_v);
    assert(u < num_v);
    return adj_mtx[v][u];
  }

  const vertex_set &neighbours(const vertex v) const {
    assert(v < num_v);
    return adj_lst[v];
  }

  vertex_set neighbourhood(const vertex v, const vertex_set &neighs) const {
    assert(v < num_v);
    vertex_set neihood;
    neihood.reserve(adj_lst[v].size());
    std::copy_if(neighs.begin(), neighs.end(), std::back_inserter(neihood),
                 [&v, this](const vertex u) { return adjacent(v, u); });
    return neihood;
  }

  graph &operator=(const graph &) = delete;
  graph &operator=(graph &&) = default;

private:
  graph(const size_type num_v, const size_type num_e)
      : num_v(num_v), num_e(num_e) {
    adj_lst.reserve(num_v);
    adj_mtx.reserve(num_v);
  }

  size_type num_v = 0, num_e = 0;

  std::vector<vertex_set> adj_lst;
  std::vector<std::vector<bool>> adj_mtx;
};
} // namespace mc

#endif // GRAPH_HPP
