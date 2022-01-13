#include "graph.hpp"

#include <iostream>
#include <algorithm>

graph::graph(const std::vector<vertex_t> &V,
			 const std::vector<edge> &E,
			 bool directed) {
  edges.reserve(V.size());
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

std::vector<vertex_t> graph::vertices() const {
  std::vector<vertex_t> verts;
  verts.reserve(edges.size());
  for (const auto &it : edges) {
	verts.emplace_back(it.first);
  }
  std::sort(verts.begin(), verts.end());
  return verts;
}

const std::vector<vertex_t> &graph::neighbors(vertex_t u) const {
  return edges.at(u);
}

std::size_t graph::degree(vertex_t u) const {
  return neighbors(u).size();
}

std::vector<std::vector<bool>> graph::adjacency_matrix() const {
  std::vector<vertex_t> verts = vertices();
  std::vector<std::vector<bool>> vec;
  vec.resize(edges.size());
  for (int i = 0; i < edges.size(); ++i) {
	vec[i].resize(edges.size());
	for (int j = 0; j < edges.size(); ++j) {
	  vertex_t u = verts.at(i), v = verts.at(j);
	  const auto &adj_lst = edges.at(u);
	  vec[i][j] = std::find(adj_lst.begin(), adj_lst.end(), v)
		!= adj_lst.end();
	}
  }
  return vec;
}

const std::unordered_map<vertex_t, std::vector<vertex_t>> &
graph::adjacency_list() const {
  return edges;
}

struct fahle_max_clique {
  struct solve {
	const std::vector<std::vector<bool>> &adj_mtx;
	std::vector<int> &sol, clique;
	solve(const std::vector<std::vector<bool>> &adj_mtx, std::vector<int> &sol)
	  : adj_mtx(adj_mtx), sol(sol) {
	  clique.reserve(adj_mtx.size());
	  std::vector<int> verts;
	  verts.reserve(adj_mtx.size());
	  for (int i = 0; i < adj_mtx.size(); ++i) {
		verts.emplace_back(i);
	  }
	  expand(std::move(verts));
	}

	std::ostream &print(int depth) {
	  return std::cout << std::string(2 * depth, ' ');
	}

	void expand(std::vector<int> verts, int depth = 0) {
	  print(depth) << " checking(" << depth << ") = {";
	  for (auto u : verts) {
		std::cout << u << ",";
	  }
	  std::cout << "}\n";

	  while (not verts.empty() && clique.size() + verts.size() > sol.size()) {
		int v = verts.back();

		print(depth) << " num verts: " << verts.size() << " ->  " << v << "\n";
		clique.emplace_back(v);

		std::vector<vertex_t> new_verts;
		for (int u : verts) {
		  print(depth)  << " --> " << v << " adjacent " << u << " = " << adj_mtx.at(v).at(u) << "\n";
		  if (adj_mtx.at(v).at(u)) {
			new_verts.emplace_back(u);
		  }
		}

		if (new_verts.empty() && clique.size() > sol.size()) {
		  sol = clique;
		  print(depth) << "solution: ";
		  for (auto u : sol) {
			print(depth) << u << " ";
		  }
		  print(depth) << "\n";
		}

		if (not new_verts.empty()) {
		  print(depth) << "expand clique of vertex: " << v << "\n";
		  expand(std::move(new_verts), depth + 1);
		}

		clique.erase(std::find(clique.begin(), clique.end(), v));
		verts.pop_back();
	  }
	  print(depth) << " checking(" << depth << ")###\n";
	  char c;
	  // std::cin >> c;
	}
  };

  fahle_max_clique(const graph *g) : g(g) {
	struct solve solver(g->adjacency_matrix(), sol);
  }

  std::vector<vertex_t> solution() const {
	std::vector<vertex_t> sol_verts, verts = g->vertices();
	std::transform(sol.begin(), sol.end(),
				   std::back_inserter(sol_verts),
				   [&verts](int i) { return verts.at(i); });
	return sol_verts;
  }


private:
  const graph *g;
  std::vector<int> sol;
};

std::vector<vertex_t> graph::fahle_max_clique() const {
  struct fahle_max_clique algo(this);
  return algo.solution();
}

bool graph::add_vertex(vertex_t v) {
  if (edges.find(v) == edges.end()) {
	edges.emplace(v);
	return true;
  }
  return false;
}

bool graph::add_edge(const edge &e) {
  add_vertex(e.from());
  add_vertex(e.to());
  auto &edges_from = edges.at(e.from());
  if (std::find(edges_from.begin(), edges_from.end(),
				e.to()) != edges_from.end()) {
	return false;
  }
  edges_from.emplace_back(e.to());
  return true;
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
	e.second.erase(std::find(e.second.begin(), e.second.end(), v));
  }
  edges.erase(it);
  return true;
}

bool graph::remove_edge(const edge &e) {
  auto it_l = edges.find(e.from());
  if (it_l == edges.end()) {
	return false;
  }
  auto it_r = std::find(it_l->second.begin(), it_l->second.end(), e.to());
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


void graph::print() const {
  for (const auto &it : edges) {
	auto size = it.second.size();
	std::cout << " vertex: " << it.first
			  << " edges: " << size << " {";
	for (const auto &e : it.second) {
	  std::cout << e << (--size ? ", " : "");
	}
	std::cout << "}\n";
  }
  std::cout << std::endl;
}

graph graph::generate() {
  std::vector<vertex_t> vertices = {1, 2, 3, 4, 5, 6, 7, 8, 9}; // vertex 0 cause a bug
  std::vector<edge> edges = {{1, 2}, {2, 3}, {3, 1}};
  return {vertices, edges, false};
}
