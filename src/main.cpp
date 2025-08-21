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

auto gen_keys(auto size) {
  std::vector<enumerator::key> keys(size);
  std::iota(keys.begin(), keys.end(), 0);
  logger::debug("V=", keys.size());
  return keys;
}

void print(const char *s, auto &key_cols) {
  std::ostringstream oss;
  oss << "using " << s << " " << key_cols.size() << " ";
  for (auto [k, col] : key_cols) {
    oss << k << "=" << col << " ";
  }
  logger::debug(oss.str());
  logger::warn(s, key_cols.highest_colour());
}

void is_valid_colouring(const char *str, auto &key_cols, enumerator &e) {
  static gc gc;

  std::pmr::unordered_map<enumerator::key, enumerator::colour> ktoc(
      gc.get_allocator());
  ktoc.reserve(key_cols.size());

  for (const auto &[u, u_col] : key_cols) {
    ktoc.emplace(u, u_col);
  }

  for (const auto &[u, u_col] : ktoc) {
    for (auto v : e.neighbours(u, {})) {
      enumerator::colour v_col = ktoc[v];
      assert(u_col != v_col, str, u_col, "should not be", v_col);
    }
  }
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

  try {
    input in;
    graph_builder builder(in);
    std::pmr::monotonic_buffer_resource vertices_pool;
    std::pmr::monotonic_buffer_resource edges_pool;
    graph g = builder.build(vertices_pool, edges_pool);
    enumerator e(g);

    auto greedy_cols = e.greedy_colour_sort(gen_keys(e.vertex_count()));
    is_valid_colouring("greedy", greedy_cols, e);
    print("greedy", greedy_cols);

    auto dsatur_cols = e.dsatur_colour_sort(gen_keys(e.vertex_count()));
    is_valid_colouring("dsatur", dsatur_cols, e);
    print("dsatur", dsatur_cols);

    return 42;

    multithreaded algo(std::move(g));

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
