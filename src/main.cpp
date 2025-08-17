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

#include "pq.hpp"

using namespace mc;

int main(int argc, char *argv[]) {
  std::set_terminate([] {
    try {
      std::rethrow_exception(std::current_exception());
    } catch (const std::runtime_error &e) {
      logger::error(e.what());
    } catch (const std::exception &e) {
      logger::error(e.what());
    } catch (...) {
      logger::error("Uncaught unknown exception!");
    }
    std::abort(); // should never be reached!
  });

  args::parse(argc, argv);

  utils::pq pq;

  pq.emplace(std::string("4"));
  assert(pq.top().str() == "4");

  pq.emplace(std::string("2"));
  assert(pq.top().str() == "2");

  pq.emplace(std::string("1"));
  assert(pq.top().str() == "1");

  pq.emplace(std::string("3"));
  assert(pq.top().str() == "1");

  pq.pop();
  assert(pq.top().str() == "2");

  pq.pop();
  assert(pq.top().str() == "3");

  pq.emplace(std::string("6"));
  assert(pq.top().str() == "3");

  pq.emplace(std::string("5"));
  assert(pq.top().str() == "3");

  pq.emplace(std::string("7"));
  assert(pq.top().str() == "3");

  pq.pop();
  assert(pq.top().str() == "4");

  pq.emplace(std::string("8"));
  assert(pq.top().str() == "4");

  pq.pop();
  assert(pq.top().str() == "5");

  pq.pop();
  assert(pq.top().str() == "6");

  pq.emplace(std::string("9"));
  assert(pq.top().str() == "6");

  pq.pop();
  assert(pq.top().str() == "7");

  pq.pop();
  assert(pq.top().str() == "8");

  pq.pop();
  assert(pq.top().str() == "9");

  pq.pop();
  assert(pq.empty());

  return 7;

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
      std::vector<graph::vertex> sorted_clique(clique.cbegin(), clique.cend());
      std::sort(sorted_clique.begin(), sorted_clique.end());
      std::ostringstream oss;
      oss << "Max Clique has " << sorted_clique.size() << " vertices { ";
      for (graph::vertex e : sorted_clique) {
        oss << e << " ";
      }
      oss << "}";
      logger::print(oss.str());
      if (args::draw) {
        algo.draw(clique);
      }
    }
  } catch (const std::exception &e) {
    std::unique_lock print_lock(logger::print_mtx);
    std::cerr << "\n\n" << e.what() << '\n';
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
