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

#include "core.hpp"
#include "mc.hpp"
#include "profiler.hpp"

using namespace std::chrono_literals;

namespace mc {
std::pmr::vector<graph::vertex> multithreaded::solve(flavour algo,
                                                     std::size_t lower_bound,
                                                     std::size_t upper_bound) {
  if (Assert(lower_bound <= upper_bound); upper_bound > 1) {
    // TODO: clean up global variables
    overall_size = lower_bound > 0 ? lower_bound - 1 : 0;
    solution(algo, upper_bound);
  } else if (Assert(G.vertex_count()); upper_bound == 1) {
    max_clique.emplace_back(0);
  }

  if (!G.is_clique(max_clique)) {
    throw std::runtime_error("NOT even a clique!!\n");
  }

  // reverse the keys back into vertices, as they were enumerated beforehand
  std::pmr::vector<graph::vertex> clique = G.to_vertex(max_clique);

  if (algo != flavour::heuristic) {
    if (max_clique.size() < lower_bound) {
      std::string err =
          "Max Clique with size=" + std::to_string(max_clique.size()) +
          " have NOT met the lower_bound=" + std::to_string(lower_bound);
      throw std::runtime_error(std::move(err));
    } else if (upper_bound == no_upper_bound &&
               (args::expect_size && max_clique.size() != args::size)) {
      std::string err =
          "Max Clique with size=" + std::to_string(max_clique.size()) +
          " is NOT maximal expected size=" + std::to_string(args::size);
      throw std::runtime_error(std::move(err));
    }
  }

  // TODO: clean up global variable assignments, since it is used recursively
  max_clique.clear();
  overall_size = 0;
  upper_bound_reached = false;
  branching_key = key();

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
  colour_sorted root = colouring_engine::colour_sort(G);

  // sort vertices based on their colours using greedy colouring
  // FIXME: re-introduce this inside the loop
  //

  //
  // inducing the neighbours out of the remaining, in fact, must be sequential

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

  thread::pool branches(args::num_threads);

  std::pmr::vector<colouring_engine> engines(memory::pool());
  engines.reserve(args::num_threads);
  for (auto i = 0u; i < args::num_threads; ++i) {
    engines.emplace_back(G);
  }

  std::string prof_name = args::filename + ".prof";

  begin = std::chrono::high_resolution_clock::now();
  profiler_start(prof_name.c_str());
  std::size_t i = root.size();
  while (!abort_search && !upper_bound_reached && i != 0) {
    i--;
    branches.exec([&, i](std::uint16_t tid) mutable {
      // make_scope_timer(branch_timer);
      if (abort_search || upper_bound_reached) {
        branches.discard_pending();
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

      auto [k, c] = root.at(i);

      auto v = G.to_vertex(k);

      if (c <= max_clique_size) {
        abort_search = true;
        branches.discard_pending();
        logger::debug(v, "has no sufficient colours. abort search.");
        return;
      }

      // make_scope_timer(sorted_neighs_timer);
      // XXX: use iterator instead to avoid blocking threads

      auto [first, last] = G.neighbours(k);
      std::pmr::vector<key> neighs(memory::pool());
      neighs.reserve(G.degree(k));

      std::pmr::vector<key> root_exclude(
          (root.get_keys().begin() + static_cast<int>(i)),
          root.get_keys().end(), memory::pool());

      std::sort(root_exclude.begin(), root_exclude.end());

      auto exclude_first = root_exclude.data();
      auto exclude_last = exclude_first + root_exclude.size();

      while (first != last && exclude_first != exclude_last) {
        if (*first < *exclude_first) {
          neighs.emplace_back(*first++);
        } else if (*exclude_first < *first) {
          exclude_first++;
        } else {
          ++first;
          ++exclude_first;
        }
      }

      std::copy(first, last, std::back_inserter(neighs));

      if (neighs.empty()) [[unlikely]] {
        abort_search = true;
        branches.discard_pending();
        logger::debug(v, "has no neighbours. abandon branching.");
        return;
      }

      colouring_engine &engine = engines[tid];
      colour_sorted sorted_neighs = engine.colour_sort(std::move(neighs));

      if (sorted_neighs.empty()) [[unlikely]] {
        abort_search = true;
        branches.discard_pending();
        logger::debug(v, "has no neighbours. abandon branching.");
        return;
      }

      // std::ostringstream oss;
      // oss << "neighs -> ";
      // for (auto u = 0ul; u < sorted_neighs.second.size(); ++u) {
      //   oss << sorted_neighs.first[u] << ":" << sorted_neighs.second[u] << "
      //   ";
      // }
      // logger::warn(oss.str());

      auto highest_colour = sorted_neighs.chromatic_num();
      if (highest_colour <= max_clique_size) {
        abort_search = true;
        branches.discard_pending();
        logger::debug("Vertex", v, "has insufficient colours. abort search.");
        return;
      }

      logger::warn("Branching on vertex", v, "with", sorted_neighs.size(),
                   "neighbours of", highest_colour, "colours");

      auto begin = std::chrono::high_resolution_clock::now();

      std::size_t num_nodes = 0;
      std::pmr::vector<key> clique(memory::pool());

      // TODO: change this to be iterative instead of recusive
      if (algo == flavour::heuristic) {
        branch_heuristic(k, k, engine, sorted_neighs, clique, max_clique_size,
                         upper_bound, num_nodes);
      } else {
        branch_exact(k, k, engine, sorted_neighs, clique, max_clique_size,
                     upper_bound, num_nodes);
      }

      // XXX: use key instead of thread_id to check for branching clique
      if (std::scoped_lock clique_lock(mtx);
          clique.size() > max_clique.size() && k == branching_key) {

        if constexpr (logger::current_level == logger::log_level::debug) {
          std::ostringstream oss;
          oss << "Found clique for " << v << " of " << clique.size()
              << " vertices { ";
          for (key u : clique) {
            oss << G.to_vertex(u) << " ";
          }
          oss << "}";
          logger::debug(oss.str());
        }

        Assert(clique.size() == overall_size);
        max_clique = std::move(clique);
      }

      std::vector<int> vff(1);
      total_branches += num_nodes;

      auto end = std::chrono::high_resolution_clock::now();
      logger::warn("Done with vertex", v, "after", num_nodes, "branches took",
                   logger::duration(begin, end));
    });
  }

  branches.join();
  profiler_stop();

  auto end = std::chrono::high_resolution_clock::now();

  logger::info(algo, "finished! found", overall_size, "vertices after",
               total_branches.load(), "branches in",
               logger::duration(begin, end));
}

bool multithreaded::enlarge_clique_size(key key, std::size_t &max_clique_size,
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

      logger::debug("Enlarged clique size to", overall_size);
    }
    max_clique_size = overall_size; // set thread local size
  }

  return found;
}

void multithreaded::branch_exact(key root, key v, colouring_engine &engine,
                                 colour_sorted &parent,
                                 std::pmr::vector<key> &clique,
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
  while (!parent.empty() && !upper_bound_reached) {
    // since the vertices were coloured "optimally", we can use them to
    // terminate a branch early when we deduce it will not lead into max
    // clique.  since the depth represents how many vertices had been
    // traversed so far, and the colour is theoretically greater than or
    // equal the potential max clique because it is a rough estimate of the
    // density of the vertex neighbourhood so as deep as we branch, the
    // colours get closer and closer to accurately indicate the size of the
    // clique we are traversing.  thus we can terminate a branch as soon as
    // we are less than the global size
    //
    // vertices we sorted in increasing degeneracy so taking the highest colour
    auto [u, c] = parent.peel();
    if (depth + c <= max_clique_size) {
      break;
    }

    // logger::debug("branching", root, "current", u, "with colour", c,
    //               "current clique", max_clique_size);
    const std::size_t prev_max_clique_size = max_clique_size;

    auto child_keys = G.neighbourhood(u, parent.get_keys());

    if (child_keys.empty() || next_depth == upper_bound) {
      if (enlarge_clique_size(root, max_clique_size, next_depth)) {
        if (next_depth == upper_bound) {
          upper_bound_reached = true; // prevent additional recursions
        }
        clique.clear();
        clique.reserve(next_depth);
        clique.emplace_back(u);
      }
    } else {
      auto child = engine.colour_sort_no_order(std::move(child_keys));
      if (!child.empty() &&
          (next_depth + child.chromatic_num()) > max_clique_size) {
        branch_exact(root, u, engine, child, clique, max_clique_size,
                     upper_bound, num_nodes, next_depth);
      }
    }

    // engine.recolour_parent(parent, c);

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
        branching_key == root) {
      clique.emplace_back(v);
    }
  }

  // logger::debug("branching off", v, "at depth", depth);
}

void multithreaded::branch_heuristic(
    key root, key v, colouring_engine &engine, colour_sorted &parent,
    std::pmr::vector<key> &clique, std::size_t &max_clique_size,
    std::size_t upper_bound, std::size_t &num_nodes, std::size_t depth) {
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

  auto [u, c] = parent.peel();
  auto child_keys = G.neighbourhood(u, parent.get_keys());

  // logger::debug("branching", root, "current", u, "with colour", c,
  //               "current clique", max_clique_size, "child_size",
  //               child_keys.size(), "parent_size", parent.size());

  if (child_keys.empty() || next_depth == upper_bound) {
    if (enlarge_clique_size(root, max_clique_size, next_depth)) {
      if (next_depth == upper_bound) {
        upper_bound_reached = true;
      }
      clique.clear();
      clique.reserve(next_depth);
      clique.emplace_back(u);
    }
  } else {
    // auto child_keys_sz = child_keys.size();
    auto child = engine.colour_sort(std::move(child_keys));
    // Assert(child.size() == child_keys_sz, "differenece in size!",
    // child.size(),
    //        child_keys_sz);
    if (!child.empty() &&
        (next_depth + child.chromatic_num()) > max_clique_size) {
      branch_heuristic(root, u, engine, child, clique, max_clique_size,
                       upper_bound, num_nodes, next_depth);
    }
  }

  // logger::debug("branching", v, "at depth", depth);

  if (prev_max_clique_size < max_clique_size && branching_key == root) {
    clique.emplace_back(v);
  }
}
} // namespace mc
