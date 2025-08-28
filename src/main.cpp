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

void pq_test() {
  utils::pq<std::string> pq;

  pq.emplace("4");
  assert(pq.top() == "4");

  pq.emplace("2");
  assert(pq.top() == "2");

  pq.emplace("1");
  assert(pq.top() == "1");

  pq.emplace("3");
  assert(pq.top() == "1");

  pq.pop();
  assert(pq.top() == "2");

  pq.pop();
  assert(pq.top() == "3");

  pq.emplace("6");
  assert(pq.top() == "3");

  pq.emplace("5");
  assert(pq.top() == "3");

  pq.emplace("7");
  assert(pq.top() == "3");

  pq.pop();
  assert(pq.top() == "4");

  pq.emplace("8");
  assert(pq.top() == "4");

  pq.pop();
  assert(pq.top() == "5");

  pq.pop();
  assert(pq.top() == "6");

  pq.emplace("9");
  assert(pq.top() == "6");

  pq.pop();
  assert(pq.top() == "7");

  pq.pop();
  assert(pq.top() == "8");

  pq.pop();
  assert(pq.top() == "9");

  pq.pop();
  assert(pq.empty());
}

void pq_destruct_test() {
  utils::pq<int> pq;

  for (int i = 0; i < 100; ++i) {
    pq.push(i);
  }
}

void exit_on_signal(int sig) { logger::error("Signal", sig, "is caught!"); }

int main(int argc, char *argv[]) {
  /* signal(SIGSEGV, exit_on_signal);
  signal(SIGFPE, exit_on_signal);
  signal(SIGILL, exit_on_signal);
  signal(SIGABRT, exit_on_signal);
  signal(SIGINT, exit_on_signal);
  signal(SIGTERM, exit_on_signal);
*/
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

  pq_test();
  pq_destruct_test();
  return 0;

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
    logger::error(e.what());
  }

  return EXIT_SUCCESS;
}
