// main.cpp
//
// Copyright (C) 2024  0x0584 (Anas)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
// USA.

#include "graph.hpp"
#include "mc.hpp"

namespace mc {
void print_clique(const std::vector<graph::vertex> &clique) {
  const std::set<graph::vertex> m(clique.cbegin(), clique.cend());
  std::ostringstream oss;
  oss << "Max Clique has " << m.size() << " vertices { ";
  for (graph::vertex e : m) {
    oss << e << " ";
  }
  oss << "}";
  log::info(oss.str());
}
} // namespace mc

int main(int argc, char *argv[]) {
  using namespace mc;

  log::setup_logger();
  args::parse(argc, argv);

  try {
    input in;

    const std::size_t vertices_buffer_size =
        in.num_v * sizeof(graph::adjacency_map::node_type);
    auto vertices_buffer = std::make_unique<std::byte[]>(vertices_buffer_size);
    std::pmr::monotonic_buffer_resource vertices_pool(vertices_buffer.get(),
                                                      vertices_buffer_size);
    const std::size_t edges_buffer_size =
        in.num_e * sizeof(graph::neighbours_set::node_type);
    auto edges_buffer = std::make_unique<std::byte[]>(edges_buffer_size);
    std::pmr::monotonic_buffer_resource edges_pool(edges_buffer.get(),
                                                   edges_buffer_size);
    graph_builder builder(in, vertices_pool, edges_pool);
    graph G = builder.build();

    multithreaded algo(std::move(G));
    for (long turn = 1; turn <= args::num_turns; ++turn) {
      log::info("Turn", turn, "/", args::num_turns);
      const std::vector<graph::vertex> clique = algo.solve(args::exec_mode);
      print_clique(clique);
      if (args::draw) {
        algo.draw(clique);
      }
    }
  } catch (const std::exception &e) {
    log::info(e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
