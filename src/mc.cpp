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

#include "mc.hpp"

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
  holder_thread_id = -1u;

  return clique;
}

void multithreaded::solution(flavour algo, std::size_t upper_bound) {
  auto begin = std::chrono::high_resolution_clock::now();

  if (algo == flavour::hybrid) {
    // TODO: handle hybrid recursion outside the function to avoid setting
    // variables on heuristic and rechecking upon them when exact branching
    if (solution(flavour::heuristic, upper_bound); upper_bound_reached) {
      auto end = std::chrono::high_resolution_clock::now();
      log::info(flavour::hybrid, "finished using", flavour::heuristic,
                "only took", log::time_diff(begin, end, log::bold));
      return;
    }

    // FIXME: prune all vertices with core number less than overall_size
  }

  std::size_t old_max_clique_size = overall_size;

  // the awaiting pool of threads has shared access to the availability vector,
  // so that if multiple threads finished simultaneously they may flag
  // termination at the same time
  std::condition_variable_any thread_pool;
  std::shared_mutex thread_mtx;
  std::vector<std::thread> threads(thread::num_threads);
  std::vector<bool> available(thread::num_threads, true);
  holder_thread_id = -1u; // TODO: handle threads in their own context
  //
  // when a vertex has no neighbours (they were pruned previously) or the
  // highest colour is less than the max clique size, then, we can terminate the
  // search since we had already sorted them in decreasing order of degeneracy
  std::atomic_bool abort_search = false;

  // vertex keys are stored in a vector, so that then would be processed
  // in-parallel, and pruned later as the algorithm proceeds
  std::vector<enumerator::key> W(E.vertex_count());
  std::iota(W.begin(), W.end(), 0);
  //
  // once a vertex is selected, it is pruned after inducing its neighbours
  const std::size_t pruned_buffer_size =
      sizeof(std::unordered_set<enumerator::key>::node_type) * E.vertex_count();
  auto pruned_buffer = std::make_unique<std::byte[]>(pruned_buffer_size);
  std::pmr::monotonic_buffer_resource pruned_pool(pruned_buffer.get(),
                                                  pruned_buffer_size);
  std::pmr::unordered_set<enumerator::key> pruned(&pruned_pool);

  //
  // inducing the neighbours out of the remaining, in fact, must be sequential
  std::mutex read_mtx;
  std::condition_variable read_flag;
  bool reading = false;

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

  auto percent_begin = std::chrono::high_resolution_clock::now();
  while (not W.empty() && not abort_search && not upper_bound_reached) {
    auto percent_end = std::chrono::high_resolution_clock::now();
    const std::size_t percent = E.vertex_count() - W.size();

#ifdef LOG
    log::info("waiting...", log::progress(percent, E.vertex_count()), "in",
              log::time_diff(begin, percent_end, log::bold));
#else // minimal console logging
    if (log::duration(percent_begin, percent_end) >= 1.) {
      log::info(log::progress(percent, E.vertex_count()), "in",
                log::time_diff(begin, percent_end, log::bold));
    }
#endif

    percent_begin = percent_end;

    {
      // block for available threads
      std::unique_lock thread_lock(thread_mtx);
      thread_pool.wait(thread_lock, [&] {
        return std::any_of(available.begin(), available.end(),
                           [](bool th) { return th; });
      });
    }

    if (abort_search || upper_bound_reached) {
      break; // search has terminated
    }

    std::size_t current_max_clique_size;
    {
      std::shared_lock clique_lock(mtx);
      current_max_clique_size = overall_size; // stamp-sync the size
    }

    if (old_max_clique_size != current_max_clique_size) {
      log::info("Max Clique of", current_max_clique_size, "vertices was found");
      old_max_clique_size = current_max_clique_size;
    }

    // sort vertices based on their colours using greedy colouring
    std::vector<enumerator::colour> colours = E.greedy_colour_sort(W);
    for (std::uint32_t thread_id = 0;
         thread_id < threads.size() && not W.empty() && not abort_search &&
         not upper_bound_reached;
         ++thread_id) {

      if (not available[thread_id]) {
        continue; // the current thread is still busy
      }

      // set the thread as busy
      if (available[thread_id].flip(); threads[thread_id].joinable()) {
        threads[thread_id].join();
      }

      if (colours.back() <= current_max_clique_size) {
#ifdef LOG
        log::print(E.key_to_vertex(W.back()),
                   "has insufficient colours. abort search!");
#endif
        abort_search = true;
        break;
      }

      enumerator::key v = W.back();
      W.pop_back();
      colours.pop_back();
      reading = true;
      threads[thread_id] = std::thread(
          [&, this](std::uint32_t thread_id, enumerator::key v,
                    std::size_t max_clique_size) {
            // since the thread can return on several conditions, it is
            // practical to use a scope destructor to ensure that threads are
            // marked as available not matter the branch
            thread::scope_dtor set_thread_available(
                [&available, &thread_id, &thread_mtx, &thread_pool] {
                  {
                    std::shared_lock thread_lock(thread_mtx);
                    available[thread_id] = true;
                  }
                  thread_pool.notify_one();
                });

            std::vector<enumerator::key> neighs = E.neighbours(v, pruned);
            pruned.emplace(v);

            {
              std::scoped_lock reading_lock(read_mtx);
              reading = false;
            }
            read_flag.notify_one();

            if (neighs.empty()) {
#ifdef LOG
              log::print(E.key_to_vertex(v),
                         "has no neighbours. abort search!");
#endif
              abort_search = true;
              return;
            }

            std::vector<enumerator::colour> colours =
                E.greedy_colour_sort(neighs);
            if (colours.back() < max_clique_size) {
#ifdef LOG
              log::print(E.key_to_vertex(v),
                         "has insufficient colours. abort search!");
#endif
              abort_search = true;
              return;
            }

#ifdef LOG
            log::print_thread(thread_id, "branching on", E.key_to_vertex(v),
                              "with", neighs.size(), "neighbours of",
                              colours.back(), "colours");
            auto begin = std::chrono::high_resolution_clock::now();
#endif

            std::size_t num_nodes = 0;
            std::vector<enumerator::key> clique;
            if (algo == flavour::heuristic) {
              branch_heuristic(thread_id, v, neighs, clique, max_clique_size,
                               upper_bound, num_nodes);
            } else {
              branch_exact(thread_id, v, neighs, colours, clique,
                           max_clique_size, upper_bound, num_nodes);
            }

            if (std::scoped_lock clique_lock(mtx);
                clique.size() > max_clique.size() &&
                thread_id == holder_thread_id) {
#ifdef LOG
              std::ostringstream oss;
              oss << "found clique for " << E.key_to_vertex(v) << " of "
                  << clique.size() << " vertices { ";
              for (enumerator::key u : clique) {
                oss << E.key_to_vertex(u) << " ";
              }
              oss << "}";
              log::print_thread(thread_id, oss.str());
#endif
              assert(clique.size() == overall_size);
              max_clique = std::move(clique);
            }

            total_branches += num_nodes;

#ifdef LOG
            auto end = std::chrono::high_resolution_clock::now();
            log::print_thread(thread_id, "done with", E.key_to_vertex(v),
                              "after", num_nodes, "branches took",
                              log::time_diff(begin, end, log::bold));
#endif
          },
          thread_id, v, current_max_clique_size);

      std::unique_lock read_lock(read_mtx);
      read_flag.wait(read_lock, [&reading]() { return not reading; });
    }
  }

#ifdef LOG
  log::info("stopped making threads");
#endif

  for (std::thread &th : threads) {
    if (th.joinable()) {
      th.join();
    }
  }
  auto end = std::chrono::high_resolution_clock::now();

#ifdef LOG
  log::info("all threads have joined");
#endif

  log::info(algo, "finished! found", overall_size, "vertices after",
            total_branches, "branches in",
            log::time_diff(begin, end, log::ansi_colours | log::bold));
}

bool multithreaded::enlarge_clique_size(std::uint32_t thread_id,
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
          holder_thread_id = thread_id;
        }
      }
      clique_shared.lock();
      found = holder_thread_id == thread_id;

#ifdef LOG
      log::print_thread(thread_id, "enlarged clique size to", overall_size);
#endif
    }
    max_clique_size = overall_size; // set thread local size
  }

  return found;
}

void multithreaded::branch_exact(std::uint32_t thread_id, enumerator::key v,
                                 std::vector<enumerator::key> &neighs,
                                 std::vector<enumerator::colour> &colours,
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
  while (not neighs.empty() && not upper_bound_reached &&
         // since the vertices were coloured "optimally", we can use them to
         // terminate a branch early when we deduce it will not lead into max
         // clique.  since the depth represents how many vertices had been
         // traversed so far, and the colour is theoretically greater than or
         // equal the potential max clique because it is a rough estimate of the
         // density of the vertex neighbourhood so as deep as we branch, the
         // colours get closer and closer to accurately indicate the size of the
         // clique we are traversing.  thus we can terminate a branch as soon as
         // we are less than the global size
         depth + colours.back() > max_clique_size) {
    const std::size_t prev_max_clique_size = max_clique_size;

    // vertices we sorted in increasing degeneracy so taking the highest colour
    enumerator::key u = neighs.back();
    neighs.pop_back();
    colours.pop_back();

    std::vector<enumerator::key> new_neighs = E.neighbourhood(u, neighs);
    if (new_neighs.empty() || next_depth == upper_bound) {
      if (enlarge_clique_size(thread_id, max_clique_size, next_depth)) {
        if (next_depth == upper_bound) {
          upper_bound_reached = true; // prevent additional recursions
        }
        clique.clear();
        clique.reserve(next_depth);
        clique.emplace_back(u);
      }
    } else {
      std::vector<enumerator::colour> new_colours =
          E.greedy_colour_sort(new_neighs);
      if (next_depth + new_colours.back() > max_clique_size) {
        branch_exact(thread_id, u, new_neighs, new_colours, clique,
                     max_clique_size, upper_bound, num_nodes, next_depth);
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
        holder_thread_id == thread_id) {
      clique.emplace_back(v);
    }
  }
}

void multithreaded::branch_heuristic(std::uint32_t thread_id, enumerator::key v,
                                     std::vector<enumerator::key> &neighs,
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
  enumerator::key u = neighs.back();
  neighs.pop_back();

  std::vector<enumerator::key> new_neighs = E.neighbourhood(u, neighs);
  if (new_neighs.empty() || next_depth == upper_bound) {
    if (enlarge_clique_size(thread_id, max_clique_size, next_depth)) {
      if (next_depth == upper_bound) {
        upper_bound_reached = true;
      }
      clique.clear();
      clique.reserve(next_depth);
      clique.emplace_back(u);
    }
  } else {
    std::vector<enumerator::colour> new_colours =
        E.greedy_colour_sort(new_neighs);
    if (next_depth + new_colours.back() > max_clique_size) {
      branch_heuristic(thread_id, u, new_neighs, clique, max_clique_size,
                       upper_bound, num_nodes, next_depth);
    }
  }

  if (prev_max_clique_size < max_clique_size && holder_thread_id == thread_id) {
    clique.emplace_back(v);
  }
}
} // namespace mc