#include <iostream>

#include "graph.hpp"

int main() {
  graph g;

  int n_verts, n_edges;
  std::cin >> n_verts >> n_edges;
  for (int i = 1; i <= n_verts; ++i) {
	g.add_vertex(i);
  }

  while (n_edges--) {
	int u, v;
	std::cin >> u >> v;
	g.add_edge({u, v});
  }

  auto &&mc = g.fahle_max_clique();
  for (vertex_t v : mc) {
	std::cout << v << " ";
  }
  std::cout << "\n";

  return 0;
}
