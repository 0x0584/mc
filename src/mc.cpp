// mc.cpp
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

#include <chrono>
#include <utility>

#include "mc.hpp"
#include "profiler.hpp"

using namespace std::chrono_literals;

namespace mc {

std::pmr::vector<graph::vertex> clique_solver::solve(flavour algo,
                                                     std::size_t lower_bound,
                                                     std::size_t upper_bound) {
  colour_sorted root = colouring_engine::colour_sort(G);

  if (Assert(lower_bound <= upper_bound); upper_bound > 1) {
    ctx.global_size = lower_bound > 0 ? lower_bound - 1 : 0;
    solution(algo, upper_bound, root);
  } else if (Assert(G.vertex_count()); upper_bound == 1) {
    ctx.global_clique.emplace_back(0);
  }

  if (!G.is_clique(ctx.global_clique)) {
    throw std::runtime_error("NOT even a clique!!\n");
  }

  std::pmr::vector<graph::vertex> clique = G.to_vertex(ctx.global_clique);

  if (algo != flavour::heuristic) {
    if (clique.size() < lower_bound) {
      std::string err =
          "Max Clique with size=" + std::to_string(clique.size()) +
          " have NOT met the lower_bound=" + std::to_string(lower_bound);
      throw std::runtime_error(std::move(err));
    } else if (upper_bound == no_upper_bound &&
               (args::expect_size && clique.size() < args::size)) {
      std::string err =
          "Max Clique with size=" + std::to_string(clique.size()) +
          " is NOT maximal expected size=" + std::to_string(args::size);
      throw std::runtime_error(std::move(err));
    }
  }

  return clique;
}

void clique_solver::solution(flavour algo, std::size_t upper_bound,
                             const colour_sorted &root) {
  auto begin = std::chrono::high_resolution_clock::now();

  if (algo == flavour::hybrid) {
    if (solution(flavour::heuristic, upper_bound, root);
        ctx.bounded.load(std::memory_order_acquire)) {
      auto end = std::chrono::high_resolution_clock::now();
      logger::info(flavour::hybrid, "finished using", flavour::heuristic,
                   "only took", logger::time_diff(begin, end, logger::bold));
      return;
    }
  }

  std::atomic_bool abort_search = false;
  thread::pool branches(args::num_threads);
  std::pmr::vector<colouring_engine> states(memory::pool());
  for (std::size_t t = 0; t < args::num_threads; ++t) {
    states.emplace_back(G);
  }

  std::ostringstream oss;
  oss << args::filename << "_" << algo << ".prof";
  std::string prof_name = oss.str();
  profiler_start(prof_name.c_str());

  begin = std::chrono::high_resolution_clock::now();
  std::atomic<std::chrono::high_resolution_clock::time_point> stamp = begin;

  std::size_t i = root.size();
  while (i != 0 && !abort_search.load(std::memory_order_acquire) &&
         !ctx.bounded.load(std::memory_order_acquire)) {
    i--;
    branches.submit([&, i](std::uint16_t tid) mutable {
      const scope_dtor abort_guard([&] {
        if (abort_search.load(std::memory_order_acquire)) [[unlikely]] {
          branches.discard_pending();
        }
      });

      auto [k, c] = root.at(i);
      if (!explore_branch(k, c)) [[unlikely]] {
        abort_search.store(true, std::memory_order_release);
        return;
      }

      std::pmr::vector<key> neighs(memory::pool());
      if (!induce_neighbours(k, root, i, neighs)) [[unlikely]] {
        abort_search.store(true, std::memory_order_release);
        return;
      }

      colouring_engine &engine = states[tid];
      auto coloured = engine.colour_sort(std::move(neighs));
      if (coloured.chromatic_num() <= ctx.snapshot_size()) [[unlikely]] {
        logger::debug("Vertex", G.to_vertex(k),
                      "has insufficient colours. abort search.");
        return;
      }

      logger::debug("Branching on vertex", G.to_vertex(k), "with",
                    coloured.size(), "neighbours of", coloured.chromatic_num(),
                    "colours with original root colour", c);

      auto stamp_val = stamp.load(std::memory_order_acquire);
      auto now = std::chrono::high_resolution_clock::now();
      if (now - stamp_val > std::chrono::seconds(5) &&
          stamp.compare_exchange_strong(stamp_val, now,
                                        std::memory_order_release,
                                        std::memory_order_relaxed)) {
        logger::info("Processed", logger::number_unit(root.size() - i),
                     "vertices so far..");
      }

      auto branch_begin = std::chrono::high_resolution_clock::now();

      if (algo == flavour::heuristic) {
        branch::runner<branch::policy::heuristic> runner(G, k, ctx, engine);
        runner.find_clique(std::move(coloured), upper_bound);
      } else {
        branch::runner<branch::policy::exact> runner(G, k, ctx, engine);
        runner.find_clique(std::move(coloured), upper_bound);
      }

      auto branch_end = std::chrono::high_resolution_clock::now();
      logger::debug("Done with vertex", G.to_vertex(k), "took",
                    logger::duration(branch_begin, branch_end));
    });
  }

  branches.join();
  profiler_stop();

  auto end = std::chrono::high_resolution_clock::now();
  logger::info(algo, "finished! found", ctx.snapshot_size(), "vertices in",
               logger::duration(begin, end));
}

} // namespace mc
