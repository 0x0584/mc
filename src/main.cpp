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

#include <climits>
#include <functional>
#include <numeric>
#include <random>
#include <set>

#include "profiler.hpp"

#include "graph.hpp"
#include "mc.hpp"
#include "thread.hpp"

using namespace mc;

int main(int argc, char *argv[]) {
  args::parse(argc, argv);

  try {
    input in;
    graph_builder builder(in);
    std::pmr::monotonic_buffer_resource vertices_pool;
    std::pmr::monotonic_buffer_resource edges_pool;
    multithreaded algo(builder.build(vertices_pool, edges_pool));
    // return 1;
    for (long turn = 1; turn <= args::num_turns; ++turn) {
      if (args::num_turns != 1)
        logger::info("Turn", turn, "/", args::num_turns);
      auto clique = algo.solve(args::exec_mode);
      auto m = std::set<graph::vertex>(clique.cbegin(), clique.cend());
      std::ostringstream oss;
      oss << "Max Clique has " << m.size() << " vertices { ";
      for (graph::vertex e : m) {
        oss << e << " ";
      }
      oss << "}";
      logger::print(oss.str());
      if (args::draw) {
        algo.draw(clique);
      }
    }
  } catch (const std::exception &e) {
    logger::info(e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
