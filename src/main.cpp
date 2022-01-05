#include "graph.hpp"

using namespace mc;

int main() {
  graph &&g = graph::generate();
  auto &&verts = g.vertices();
  g.print();
  for (auto u : verts) {
	for (auto v : verts) {
	  if (g.remove_edge({u, v}))
		g.print();
	}
	std::cout << " --- \n";
  }
  return 0;
}
