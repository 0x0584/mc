#include "graph.hpp"

namespace mc {
  graph::graph(const std::vector<vertex_t> &V,
			   const std::vector<edge> &E,
			   bool directed) {
	for (const auto &v : V) {
	  add_vertex(v);
	}
	for (const auto &e : E) {
	  if (directed) {
		add_edge(e);
	  } else {
		add_edge_undirected(e);
	  }
	}
  }

  bool graph::add_vertex(vertex_t v) {
	if (edges.find(v) == edges.end()) {
	  edges.emplace(v);
	  return true;
	}
	return false;
  }

  bool graph::add_edge(const edge &e) {
	add_vertex(e.to());
	return edges[e.from()].emplace(e.to(), e.from()).second;
  }

  void graph::add_edge_undirected(const edge &e) {
	add_edge(e);
	add_edge(e.residual());
  }

  bool graph::remove_vertex(vertex_t v) {
	auto it = edges.find(v);
	if (it == edges.end()) {
	  return false;
	}
	for (auto &e : edges) {
	  e.second.erase(v);
	}
	edges.erase(it);
	return true;
  }

  bool graph::remove_edge(const edge &e) {
	auto it_l = edges.find(e.from());
	if (it_l == edges.end()) {
	  return false;
	}
	auto it_r = it_l->second.find(e.to());
	if (it_r == it_l->second.end()) {
	  return false;
	}
	it_l->second.erase(it_r);
	return true;
  }

  void graph::remove_edge_undirected(const edge &e) {
	remove_edge(e);
	remove_edge(e.residual());
  }

  std::vector<vertex_t> graph::vertices() const {
	std::vector<vertex_t> verts;
	verts.reserve(edges.size());
	for (const auto &it : edges) {
	  verts.emplace_back(it.first);
	}
	return verts;
  }

  // std::vector<edge> graph::edges() const {

  // }


  void graph::print() const {
	for (const auto &it : edges) {
	  auto size = it.second.size();
	  std::cout << " vertex_t: " << it.first
				<< " edges: " << size << " {";
	  for (const auto &e : it.second) {
		std::cout << e.first << (--size ? ", " : "");
	  }
	  std::cout << "}\n";
	}
	std::cout << std::endl;
  }

  graph graph::generate() {
	std::vector<vertex_t> vertices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	std::vector<edge> edges = {{1, 2}, {1, 3}, {2, 3}, {3, 5},
							   {4, 1}, {4, 2}, {5, 6}, {7, 1},
							   {8, 1}, {9, 1}, {9, 2}, {9, 3}};
	return {vertices, edges, true};
  }
}
