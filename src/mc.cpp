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

#include <gperftools/profiler.h>

#include <barrier>
#include <chrono>
#include <condition_variable>
#include <memory_resource>
#include <mutex>

#include <thread>
#include <utility>

#include <algorithm>
#include <numeric>
#include <unordered_set>

#include "mc.hpp"

using namespace std::chrono_literals;

namespace mc {
std::vector<graph::vertex> multithreaded::solve(flavour algo,
                                                std::size_t lower_bound,
                                                std::size_t upper_bound) {
  if (assert(lower_bound <= upper_bound); upper_bound > 1) {
    // TODO: clean up global variables
    overall_size = lower_bound > 0 ? lower_bound - 1 : 0;
    solution(algo, upper_bound);
  } else if (assert(E.vertex_count()); upper_bound == 1) {
    max_clique.emplace_back(0);
  }

  // reverse the keys back into vertices, as they were enumerated beforehand
  std::vector<graph::vertex> clique = E.unfold_keys(max_clique);

  if (not E.is_clique(max_clique)) {
    throw std::runtime_error("NOT even a clique!!\n");
  }

  if (algo != flavour::heuristic) {
    if (max_clique.size() < lower_bound) {
      const std::string err =
          " Max Clique with size=" + std::to_string(max_clique.size()) +
          " have NOT met the lower_bound=" + std::to_string(lower_bound) + "\n";
      throw std::runtime_error(err.c_str());
    } else if (upper_bound == no_upper_bound &&
               (args::expect_size && max_clique.size() != args::size)) {
      const std::string err =
          " Max Clique with size=" + std::to_string(max_clique.size()) +
          " is NOT maximal expected size=" + std::to_string(args::size) + "\n";
      throw std::runtime_error(err.c_str());
    }
  }

  // TODO: clean up global variable assignments, since it is used recursively
  max_clique.clear();
  overall_size = 0;
  upper_bound_reached = false;
  branching_key = enumerator::null_key;

  return clique;
}

void multithreaded::solution(flavour algo, std::size_t upper_bound) {
  auto begin = std::chrono::high_resolution_clock::now();

  if (algo == flavour::hybrid) {
    // TODO: handle hybrid recursion outside the function to avoid setting
    // variables on heuristic and rechecking upon them when exact branching
    if (solution(flavour::heuristic, upper_bound); upper_bound_reached) {
      auto end = std::chrono::high_resolution_clock::now();
      logger::info(flavour::hybrid, "finished using", flavour::heuristic,
                   "only took", logger::time_diff(begin, end, logger::bold));
      return;
    }

    // FIXME: prune all vertices with core number less than overall_size
  }

  std::size_t old_max_clique_size = overall_size;

  // when a vertex has no neighbours (they were pruned previously) or the
  // highest colour is less than the max clique size, then, we can terminate the
  // search since we had already sorted them in decreasing order of degeneracy
  std::atomic_bool abort_search = false;

  // vertex keys are stored in a vector, so that then would be processed
  // in-parallel, and pruned later as the algorithm proceeds
  enumerator::sorted_keys sorted_keys;
  {
    std::vector<enumerator::key> keys(E.vertex_count());
    std::iota(keys.begin(), keys.end(), 0);
    sorted_keys = E.greedy_colour_sort(std::move(keys));
  }
  // sort vertices based on their colours using greedy colouring
  // FIXME: re-introduce this inside the loop
  //
  std::pmr::unordered_set<enumerator::key> pruned;
  //
  // inducing the neighbours out of the remaining, in fact, must be sequential
  std::mutex pruned_mtx;

  // while branches run in-parallel, each might not visit all vertices
  // thus we define a branching as follows:
  //
  //   1. select a vertex that has not been pruned
  //   2. induce its neighbourhood of the remaining vertices
  //   3. colour the neighbours
  //   4. recurs on all vertices (until we have an empty set of neighbours) only
  //      if the maximum colour is greater than the depth of the branch
  //   5. if the depth is greater than the current max clique size, backtrack
  //      and gather vertices in the current max clique
  //
  // since the algorithm is recursive, each callback is a branching
  std::atomic_size_t total_branches = 0;
  E.reset_cache_hits();

#ifndef NDEBUG
  ProfilerStart("max-clique.prof");
#endif

  thread::pool branches(args::num_threads);

  begin = std::chrono::high_resolution_clock::now();

  std::barrier prune_barrier(2);
  while (not sorted_keys.empty() && not abort_search &&
         not upper_bound_reached) {
    auto percent_begin = std::chrono::high_resolution_clock::now();
    auto percent = E.vertex_count() - sorted_keys.size();

    branches.exec([this, &abort_search, &old_max_clique_size, &pruned, &algo,
                   &upper_bound, &total_branches, &prune_barrier,
                   key_colour =
                       sorted_keys.pop()](std::uint16_t task_id) mutable {
      bool arrived = false;
      scope_dtor barrier([&arrived, &prune_barrier] {
        if (not arrived) {
          prune_barrier.arrive_and_wait();
        }
      });

      if (upper_bound_reached || abort_search) {
        return;
      }

      std::size_t max_clique_size;
      {
        std::shared_lock clique_lock(mtx);
        max_clique_size = overall_size;
      }

      if (old_max_clique_size != max_clique_size) {
        old_max_clique_size = max_clique_size;
        logger::debug("New Clique with", max_clique_size,
                      "vertices was found!");
      }

      const auto &[key, colour] = key_colour;
      auto v = E.key_to_vertex(key);

      if (colour <= max_clique_size) {
        logger::debug(v, "has no sufficient colours. abort search.");
        abort_search = true;
        return;
      }

      auto neighs = E.neighbours(key, pruned);
      pruned.emplace(key);
      prune_barrier.arrive_and_wait();
      arrived = true;

      if (neighs.empty()) {
        logger::debug(v, "has no neighbours. abandon branching.");
        return;
      }

      auto sorted_neighs = E.greedy_colour_sort(std::move(neighs));
      if (sorted_neighs.highest_colour() < max_clique_size) {
        abort_search = true;
        logger::debug("Vertex", v, "has insufficient colours. abort search.");
        return;
      }

      logger::debug("Branching on vertex", v, "with", sorted_neighs.size(),
                    "neighbours of", sorted_neighs.highest_colour(), "colours");

      auto begin = std::chrono::high_resolution_clock::now();

      std::size_t num_nodes = 0;
      std::vector<enumerator::key> clique;

      if (algo == flavour::heuristic) {
        branch_heuristic(key, key, sorted_neighs, clique, max_clique_size,
                         upper_bound, num_nodes);
      } else {
        branch_exact(key, key, sorted_neighs, clique, max_clique_size,
                     upper_bound, num_nodes);
      }

      // FIXME: use key instead of thread_id to check for branching clique
      if (std::scoped_lock clique_lock(mtx);
          clique.size() > max_clique.size() && key == branching_key) {

        if constexpr (logger::current_level == logger::log_level::debug) {
          std::ostringstream oss;
          oss << "Found clique for " << v << " of " << clique.size()
              << " vertices { ";
          for (enumerator::key u : clique) {
            oss << E.key_to_vertex(u) << " ";
          }
          oss << "}";
          logger::debug(oss.str());
        }

        assert(clique.size() == overall_size);
        max_clique = std::move(clique);
      }

      total_branches += num_nodes;

      auto end = std::chrono::high_resolution_clock::now();
      logger::debug("Done with vertex", v, "after", num_nodes, "branches took",
                    logger::time_diff(begin, end, logger::bold), "with",
                    E.get_cache_hits(), "cache hits");
    });

    prune_barrier.arrive_and_wait();

    auto percent_end = std::chrono::high_resolution_clock::now();

    if constexpr (logger::current_level == logger::log_level::debug) {
      if (logger::duration(percent_begin, percent_end) >= 1.) {
        logger::debug(logger::progress(percent, E.vertex_count()), "in",
                      logger::time_diff(begin, percent_end, logger::bold));
      }
    }
  }

  branches.join();

  auto end = std::chrono::high_resolution_clock::now();

#ifndef NDEBUG
  ProfilerStop();
#endif

  logger::info(
      algo, "finished! found", overall_size, "vertices after", total_branches,
      "branches and", E.get_cache_hits(), "cache hits with ratio of",
      (double(E.get_cache_hits()) / total_branches), "in",
      logger::time_diff(begin, end, logger::ansi_colours | logger::bold));
}

bool multithreaded::enlarge_clique_size(enumerator::key key,
                                        std::size_t &max_clique_size,
                                        std::size_t depth) {
  bool found = false;
  {
    std::shared_lock clique_shared(mtx);
    if (depth > overall_size) {
      clique_shared.unlock();
      {
        // ensure sequential (despite it being out of order) execution
        std::scoped_lock clique_lock(mtx);
        // only after that we save the old value for further comparison
        const std::size_t old_overall_size = overall_size;
        // ensuring that we take the correct (intended) maximum value of both
        // so if indeed we found a great size, then we set this thread as the
        // holder of the current maximum clique
        if (overall_size = std::max(depth, overall_size);
            old_overall_size != overall_size) {
          branching_key = key;
        }
      }
      clique_shared.lock();
      found = branching_key == key;

#ifdef LOG
      logger::print_thread(thread_id, "enlarged clique size to", overall_size);
#endif
    }
    max_clique_size = overall_size; // set thread local size
  }

  return found;
}

void multithreaded::branch_exact(enumerator::key key, enumerator::key v,
                                 enumerator::sorted_keys &sorted_neighs,
                                 std::vector<enumerator::key> &clique,
                                 std::size_t &max_clique_size,
                                 std::size_t upper_bound,
                                 std::size_t &num_nodes, std::size_t depth) {
  num_nodes++;

  {
    // updating the local max clique often to ensure pruning vertices as
    // accurately as possible in the context of all the running threads
    std::shared_lock clique_lock(mtx);
    max_clique_size = overall_size;
  }

  const std::size_t next_depth = depth + 1;
  //  const std::size_t neighs_before = sorted_neighs.size();
  while (not sorted_neighs.empty() && not upper_bound_reached &&
         // since the vertices were coloured "optimally", we can use them to
         // terminate a branch early when we deduce it will not lead into max
         // clique.  since the depth represents how many vertices had been
         // traversed so far, and the colour is theoretically greater than or
         // equal the potential max clique because it is a rough estimate of the
         // density of the vertex neighbourhood so as deep as we branch, the
         // colours get closer and closer to accurately indicate the size of the
         // clique we are traversing.  thus we can terminate a branch as soon as
         // we are less than the global size
         depth + sorted_neighs.highest_colour() > max_clique_size) {
    const std::size_t prev_max_clique_size = max_clique_size;

    // vertices we sorted in increasing degeneracy so taking the highest colour
    enumerator::key u = sorted_neighs.pop_key_with_highest_colour();

    std::vector<enumerator::key> new_neighs =
        E.neighbourhood(u, sorted_neighs.keys());
    if (new_neighs.empty() || next_depth == upper_bound) {
      if (enlarge_clique_size(key, max_clique_size, next_depth)) {
        if (next_depth == upper_bound) {
          upper_bound_reached = true; // prevent additional recursions
        }
        clique.clear();
        clique.reserve(next_depth);
        clique.emplace_back(u);
      }
    } else {
      auto new_sorted_neighs = E.greedy_colour_sort(std::move(new_neighs));
      if (next_depth + new_sorted_neighs.highest_colour() > max_clique_size) {
        branch_exact(key, u, new_sorted_neighs, clique, max_clique_size,
                     upper_bound, num_nodes, next_depth);
      }
    }

    // when we reach a leaf, when recursion terminates, we had already saved the
    // old size so we can determine if a clique had been found for the current
    // vertex v since the whole recursion only concerns a single vertex
    if (prev_max_clique_size < max_clique_size &&
        // although on a surface level, access without locking the mutex seems
        // like a fatal data-race, but in fact it is correct in this specific
        // context! since we only perform a read-op and only the holder thread
        // is guaranteed to have the correct value: which is exactly what we are
        // trying to achieve.  practically (as far as I had measured) locking
        // the mutex would x2 the runtime with not clear benefit
        branching_key == key) {
      clique.emplace_back(v);
    }
  }

  // logger::print_thread(thread_id, "branching off", E.key_to_vertex(v),
  // "with",
  //                   neighs_before - neighs.size(),
  //                   "remaining neighbours with max_colour", colours.back(),
  //                   "colours and depth", depth);
}

void multithreaded::branch_heuristic(enumerator::key key, enumerator::key v,
                                     enumerator::sorted_keys &sorted_neighs,
                                     std::vector<enumerator::key> &clique,
                                     std::size_t &max_clique_size,
                                     ::size_t upper_bound,
                                     std::size_t &num_nodes,
                                     std::size_t depth) {
  if (upper_bound_reached) {
    return;
  }

  num_nodes++;

  {
    std::shared_lock clique_lock(mtx);
    max_clique_size = overall_size;
  }

  const std::size_t prev_max_clique_size = max_clique_size;
  const std::size_t next_depth = depth + 1;

  // the essence of the heuristic is instead of traversing all neighbours, we
  // only pick the most promising one: the vertex with the highest colour
  enumerator::key u = sorted_neighs.pop_key_with_highest_colour();

  std::vector<enumerator::key> new_neighs =
      E.neighbourhood(u, sorted_neighs.keys());
  if (new_neighs.empty() || next_depth == upper_bound) {
    if (enlarge_clique_size(key, max_clique_size, next_depth)) {
      if (next_depth == upper_bound) {
        upper_bound_reached = true;
      }
      clique.clear();
      clique.reserve(next_depth);
      clique.emplace_back(u);
    }
  } else {
    auto new_sorted_neighs = E.greedy_colour_sort(std::move(new_neighs));
    if (next_depth + new_sorted_neighs.highest_colour() > max_clique_size) {
      branch_heuristic(key, u, new_sorted_neighs, clique, max_clique_size,
                       upper_bound, num_nodes, next_depth);
    }
  }

  if (prev_max_clique_size < max_clique_size && branching_key == key) {
    clique.emplace_back(v);
  }
}
} // namespace mc
