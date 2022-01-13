#pragma once

#include <unordered_map>
#include <vector>

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
  graph();
  graph(const graph &) = delete;
  graph(const std::vector<vertex_t> &V,
		const std::vector<edge> &E,
		bool directed = true);

  std::vector<vertex_t> vertices() const;
  const std::vector<vertex_t> &neighbors(vertex_t u) const;
  std::size_t num_vertices() const { return edges.size(); }
  std::size_t degree(vertex_t u) const;
  std::vector<std::vector<bool>> adjacency_matrix() const;
  const std::unordered_map<vertex_t, std::vector<vertex_t>> &adjacency_list() const;

  std::vector<vertex_t> fahle_max_clique() const;

  bool add_vertex(vertex_t v);
  bool add_edge(const edge &e);
  void add_edge_undirected(const edge &e);

  bool remove_vertex(vertex_t v);
  bool remove_edge(const edge &e);
  void remove_edge_undirected(const edge &e);

  void print() const;

  graph &operator=(const graph &) = delete;

  static graph generate();

private:
  std::unordered_map<vertex_t, std::vector<vertex_t>> edges;
};
