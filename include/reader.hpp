#ifndef READER_HPP
#define READER_HPP

#include "graph.hpp"

namespace mc {
namespace order_vertices {
struct order_by_neighbours;
struct order_by_neighbourhood;
struct order_by_degeneracy;

using vertex_order = std::variant<order_by_neighbours, order_by_neighbourhood,
                                  order_by_degeneracy>;
} // namespace order_vertices

struct reader {
  using vertex = size_type;
  using vertex_order = order_vertices::vertex_order;
  using vertex_set = std::vector<vertex>;
  using vertex_map = std::unordered_map<vertex, graph::vertex>;

  reader() = delete;
  reader(const reader &) = delete;

  reader(const vertex_order &order, std::istream &in = std::cin);

  size_type num_vertices() const { return g.num_v; }
  size_type num_edges() const { return g.num_e; }

  const graph &get() const { return g; }

  graph::vertex vertex_to_key(const vertex &v) const {
    assert(v2k.count(v));
    return v2k.at(v);
  }

  const vertex &key_to_vertex(const graph::vertex v) const {
    assert(v < g.num_v);
    return k2v[v];
  }

  reader &operator=(const reader &) = delete;

  friend std::ostream &operator<<(std::ostream &os, const reader &gr);

private:
  graph g;

  vertex_set k2v;
  vertex_map v2k;
};

using adjacency_list =
    std::unordered_map<reader::vertex, std::set<reader::vertex>>;
using adjacency_list_sorted =
    std::vector<std::pair<reader::vertex, std::set<reader::vertex>>>;

namespace order_vertices {
struct order_by_neighbours {
  adjacency_list_sorted operator()(adjacency_list &adj_lst_orig) const;
};

struct order_by_neighbourhood {
  adjacency_list_sorted operator()(adjacency_list &adj_lst_orig) const;
};

struct order_by_degeneracy {
  adjacency_list_sorted operator()(adjacency_list &adj_lst_orig) const;
};

struct dispatcher {
  explicit dispatcher(adjacency_list &adj_lst_orig)
      : adj_lst_orig(adj_lst_orig) {}

  dispatcher() = delete;
  dispatcher(const dispatcher &) = delete;
  dispatcher(dispatcher &&) = delete;

  template <typename Order> adjacency_list_sorted operator()(Order &&sort) {
    MAKE_TIMER(dispatcher);
    return sort(adj_lst_orig);
  }

  dispatcher &operator=(const dispatcher &) = delete;
  dispatcher &operator=(dispatcher &&) = delete;

private:
  adjacency_list &adj_lst_orig;
};

extern order_by_neighbours by_neighbours;
extern order_by_neighbourhood by_neighbourhood;
extern order_by_degeneracy by_degeneracy;
} // namespace order_vertices
} // namespace mc
#endif // READER_HPP
