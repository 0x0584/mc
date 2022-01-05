#pragma once

#include <iostream>

#include <unordered_map>
#include <vector>

namespace mc {
  using vertex_t = int;

  struct edge {
	using map = std::unordered_map<vertex_t, vertex_t>;

	edge &operator=(const edge &) = delete;

	edge(vertex_t u, vertex_t v) : s(u), t(v) {}
	edge(const map::value_type &e) : s(e.first), t(e.second) {}
	edge(const std::pair<vertex_t, vertex_t> &e)
	  : s(e.first), t(e.second) {}

	edge residual() const { return {t, s}; }

	vertex_t from() const { return s; }
	vertex_t to() const { return t; }

  private:
	vertex_t s, t;
  };

  class graph {
	std::unordered_map<vertex_t, edge::map> edges;

  public:

	graph();
	graph(const std::vector<vertex_t> &V,
		  const std::vector<edge> &E,
		  bool directed = true);

	graph(const graph &) = delete;

	bool add_vertex(vertex_t v);
	bool add_edge(const edge &e);
	void add_edge_undirected(const edge &e);

	bool remove_vertex(vertex_t v);
	bool remove_edge(const edge &e);
	void remove_edge_undirected(const edge &e);

	std::vector<vertex_t> vertices() const;
	//	std::vector<edge> edges() const;

	void print() const;

	graph &operator=(const graph &) = delete;

	static graph generate();
  };
}
