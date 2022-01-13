#include "graph.hpp"

#include <iostream>

void print(const graph &g) {
  std::cout << "Adjacency List\n\n";
  g.print();
  std::cout << "Adjacency Matrix\n\n  ";
  auto &&verts = g.vertices();
  for (auto v : verts) {
	std::cout << v << " ";
  }
  std::cout << "\n";
  auto &&adj_mtx = g.adjacency_matrix();
  for (vertex_t u = 0; u < verts.size(); ++u) {
	std::cout << verts[u] << " ";
	for (auto col : adj_mtx.at(u)) {
	  std::cout << col << " ";
	}
	std::cout << "\n";
  }
  std::cout << "\n";
}

int main() {
  graph &&g = graph::generate();
  auto &&verts = g.vertices();
  print(g);
  auto &&mc = g.fahle_max_clique();
  std::cout << " ---------------- \n";
  for (auto u : mc) {
	std::cout << u << " ";
  }
  std::cout << "\n ---------------- \n";

  // print(g);
  // for (auto u : verts) {
  // 	for (auto v : verts) {
  // 	  if (g.remove_edge({u, v})) {
  // 		print(g);
  // 	  }
  // 	}
  // 	std::cout << " --- \n";
  // }


  return 0;
}
