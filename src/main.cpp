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

#include "colour.hpp"
#include "graph.hpp"
#include "mc.hpp"

#include <sys/stat.h>
#include <unistd.h>

static inline bool is_dev_null(int fd) {
  struct stat fd_stat;
  struct stat dev_null_stat;

  // Get status of the file descriptor
  if (fstat(fd, &fd_stat) == -1) {
    perror("fstat");
    return false;
  }

  // Get status of /dev/null
  if (stat("/dev/null", &dev_null_stat) == -1) {
    perror("stat");
    return false;
  }

  // Compare device and inode numbers
  return (fd_stat.st_dev == dev_null_stat.st_dev) &&
         (fd_stat.st_ino == dev_null_stat.st_ino);
}

const bool no_stdout = is_dev_null(1);
const bool no_stderr = is_dev_null(2);

using namespace mc;

std::string colour_sorted_str(const colour_sorted &cs) {
  std::ostringstream oss;
  oss << "(";
  if (!cs.empty()) {
    for (auto i = 0ul; i < cs.size(); ++i) {
      oss << " " << cs.key_at(i) << "=" << cs.colour_at(i);
    }
  }
  oss << " )";
  return oss.str();
}

std::string vector_str(const std::pmr::vector<key> &v) {
  std::ostringstream oss;
  oss << "(";
  for (auto &u : v) {
    oss << " " << u;
  }
  oss << " )";
  return oss.str();
}

std::pmr::vector<key> neighbours(const graph &g, key u) {
  auto [first, last] = g.neighbours(u);
  std::pmr::vector<key> v{first, last};
  // logger::print(u, "orig neighbours", vector_str(v));
  return v;
}

void simulate_branching(const graph &g, colouring_engine &engine, key u,
                        colour_sorted &parent, std::size_t depth) {
  std::string depth_str(depth, ' ');

  for (int i = 1; !parent.empty(); ++i) {
    logger::info(depth_str, "depth=", depth, "parent=", u,
                 "chromatic=", parent.chromatic_num(),
                 "coloured=", colour_sorted_str(parent));

    auto [v, c] = parent.peel();
    auto R = g.neighbourhood(v, parent.get_keys());

    auto R_str = vector_str(R);

    colour_sorted child = engine.colour_sort(std::move(R));

    logger::debug(depth_str, "i=", i, u, "is parent (peeled=", v,
                  "chromatic=", c, ") child=", v, "neighbourhood=", R_str,
                  "child=", colour_sorted_str(child));

    simulate_branching(g, engine, v, child, depth + 1);
  }

  logger::info(depth_str, "done depth=", depth, "parent=", u);
}

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

  input in;
  graph_builder builder(in);
  graph g = builder.build(args::undirected);

  // key u = 20;

  // auto u_neighs = neighbours(g, u);
  // std::sort(u_neighs.begin(), u_neighs.end());
  // logger::print(u, "neighbours", vector_str(u_neighs));

  // colouring_engine engine(g);

  // auto cs = engine.colour_sort(std::move(u_neighs));

  // do {
  //   u_neighs = cs.get_keys();
  //   cs = engine.colour_sort(std::move(u_neighs));
  //   logger::print(u, "colours", colour_sorted_str(cs));
  //   cs.peel();
  // } while (!cs.empty());

  // return 42;

  /*
  key u = 20;
  colouring_engine engine(g);
  colour_sorted u_coloured = engine.colour_sort_first(neighbours(g, u));
  is_valid_colouring(g, u_coloured);

  logger::info("ROOT parent=", u, " coloured=", colour_sorted_str(u_coloured));
  simulate_branching(g, engine, u, u_coloured, 1);

  return 42;
*/
  // auto keys = gen_keys(g.vertex_count());

  // auto t0 = std::chrono::high_resolution_clock::now();
  // colour_sorted coloured_g = colouring_engine::colour_sort_graph(g);
  // auto t1 = std::chrono::high_resolution_clock::now();
  // logger::warn("colouring_took", logger::duration(t0, t1));
  // if (!no_stdout) {
  //   std::ostringstream oss;
  //   oss << "chromatic number: " << coloured_g.chromatic_num() << "\n";
  //   for (auto i = 0ul; i < coloured_g.size(); ++i) {
  //     oss << g.to_vertex(coloured_g.key_at(i));
  //     oss << "(" << coloured_g.colour_at(i) << ") ";
  //   }
  //   logger::print(oss.str());
  // }
  // is_valid_colouring(g, coloured_g);
  // return 42;

  clique_solver algo(g);
  for (long turn = 1; turn <= args::num_turns; ++turn) {
    if (args::num_turns != 1) {
      logger::info("Turn", turn, "/", args::num_turns);
    }

    auto clique = algo.solve(args::exec_mode);
    std::pmr::vector<graph::vertex> sorted_clique(
        clique.cbegin(), clique.cend(), memory::pool());
    std::sort(sorted_clique.begin(), sorted_clique.end());
    std::ostringstream oss;
    oss << "Max Clique has " << sorted_clique.size() << " vertices { ";
    for (const graph::vertex &v : sorted_clique) {
      oss << v << " ";
    }
    oss << "}";
    logger::print(oss.str());

    if (args::draw) {
      // algo.draw(clique);
    }
  }

  return EXIT_SUCCESS;
}
