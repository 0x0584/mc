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

  // g.add_edge_undirected({1, 2});
  // g.add_edge_undirected({1, 4});
  // g.add_edge_undirected({1, 3});
  // g.add_edge_undirected({1, 8});

  // g.add_edge_undirected({2, 6});
  // g.add_edge_undirected({2, 5});
  // g.add_edge_undirected({2, 4});
  // g.add_edge_undirected({2, 7});

  // g.add_edge_undirected({3, 4});
  // g.add_edge_undirected({3, 4});
  // g.add_edge_undirected({3, 6});
  // g.add_edge_undirected({3, 5});
  // g.add_edge_undirected({3, 10});

  // g.add_edge_undirected({4, 6});
  // g.add_edge_undirected({4, 7});
  // g.add_edge_undirected({4, 5});
  // g.add_edge_undirected({4, 9});
  // g.add_edge_undirected({4, 10});

  // g.add_edge_undirected({6, 8});

  // g.add_edge_undirected({5, 6});
  // g.add_edge_undirected({5, 7});
  // g.add_edge_undirected({5, 8});

  // g.add_edge_undirected({7, 9});
  // g.add_edge_undirected({7, 9});
  // g.add_edge_undirected({7, 10});

  // g.add_edge_undirected({8, 9});

  // g.add_edge_undirected({9, 10});

  auto &&mc = g.fahle_max_clique();
  std::sort(mc.begin(), mc.end());
  std::cout << "SIZE = " << mc.size() << " { ";
  for (vertex_t v : mc) {
    std::cout << v << " ";
  }
  std::cout << "}\n";

  return 0;
}
