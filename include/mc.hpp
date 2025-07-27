// mc.hpp
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

#ifndef MAXCLIQUE_HPP
#define MAXCLIQUE_HPP

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <execution>
#include <functional>
#include <future>
#include <memory_resource>
#include <mutex>
#include <shared_mutex>
#include <thread>
#include <utility>

using namespace std::chrono_literals;

#include <algorithm>
#include <deque>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// #define NDEBUG

#include "enumerator.hpp"
#include "thread.hpp"

namespace mc {
class multithreaded {
public:
  static inline const std::size_t maximum_bound = -1u;

  enumerator E;

private:
  // only a single mutex is used to handle the max_clique and its global size,
  // in addition to which thread because the max_clique is updated only when we
  // finished branching, and the size is updated during the branching.  it is
  // also shared since most of the time we just want a read-op, so it is optimal
  // to use std::shared_mtx
  std::shared_mutex mtx;
  //
  // the size of the largest clique found so far across all the running threads,
  // however, it is updated separately from the actual max_clique with the depth
  // of the branch rather than actually counting the clique vertices
  std::size_t overall_size = 0;
  //
  // after branch termination, if the current thread had found the largest one
  // so far amongst all the running threads (even if they are still running
  std::vector<enumerator::key> max_clique;
  //
  // hence, we can set the max clique few times and avoid unnecessary
  //  assignments of cliques from several threads, at least in most cases
  std::uint32_t holder_thread_id = -1u;

  // terminate the algorithm early if the depth matches the bound
  std::atomic_bool upper_bound_reached = false;

  void solution(flavour algo, std::size_t upper_bound);

  bool enlarge_clique_size(std::uint32_t thread_id,
                           std::size_t &max_clique_size, std::size_t depth);

  void branch_exact(std::uint32_t thread_id, enumerator::key v,
                    enumerator::sorted_keys &sorted_neighs,
                    std::vector<enumerator::key> &clique,
                    std::size_t &max_clique_size, std::size_t upper_bound,
                    std::size_t &num_nodes, std::size_t depth = 1);

  void branch_heuristic(std::uint32_t thread_id, enumerator::key v,
                        enumerator::sorted_keys &sorted_neighs,
                        std::vector<enumerator::key> &clique,
                        std::size_t &max_clique_size, std::size_t upper_bound,
                        std::size_t &num_nodes, std::size_t depth = 1);

public:
  static inline std::size_t no_upper_bound = -1u;

  explicit multithreaded(graph G) : E(std::move(G)) {
    log::info("Number of available Threads", thread::num_threads);
  }

  ~multithreaded() { log::info("~multithreaded()"); }

  std::vector<graph::vertex>
  solve(flavour algo = flavour::exact,
        // the expected behaviour is (as far as I have tested) the function call
        // with be launched with the up-to-date values, even though the it seems
        // to be at compile time, it is dynamic initialisation
        std::size_t lower_bound = args::lower_bound,
        std::size_t upper_bound = args::upper_bound);

  void draw(const std::vector<graph::vertex> &clique) { E.draw(clique); }
};
} // namespace mc
#endif // MAXCLIQUE_HPP
