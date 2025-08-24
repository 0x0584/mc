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
#include <shared_mutex>
#include <vector>

// #define NDEBUG

#include "cache.hpp"
#include "graph.hpp"

namespace mc {
class multithreaded {
  /*struct branch {
    branch(const enumerator &E,
           const std::pmr::unordered_set<enumerator::key> &pruned)
        : E(E), pruned(pruned) {}

  protected:
    std::size_t num_nodes = 0;
    std::pmr::vector<enumerator::key> clique;

    const enumerator &E;
    const std::pmr::unordered_set<enumerator::key> &pruned;
  };

  struct branch_heuristic : branch {};

  struct branch_exact : branch {};*/

public:
  static inline const std::size_t maximum_bound = -1u;

private:
  const graph &G;

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
  std::pmr::vector<graph::key> max_clique;
  //
  // hence, we can set the max clique few times and avoid unnecessary
  //  assignments of cliques from several threads, at least in most cases
  graph::key branching_key;

  // terminate the algorithm early if the depth matches the bound
  std::atomic_bool upper_bound_reached = false;

  void solution(flavour algo, std::size_t upper_bound);

  bool enlarge_clique_size(graph::key key, std::size_t &max_clique_size,
                           std::size_t depth);

  void branch_exact(
      graph::key key, graph::key v,
      std::pair<std::pmr::vector<graph::key>, std::pmr::vector<colour>>
          &sorted_neighs,
      std::pmr::vector<graph::key> &clique, std::size_t &max_clique_size,
      std::size_t upper_bound, std::size_t &num_nodes,
      lru_cache<std::pmr::vector<graph::key>, std::pmr::vector<colour>> &cache,
      std::size_t &cache_hits, std::size_t depth = 1);

  void branch_heuristic(
      graph::key key, graph::key v,
      std::pair<std::pmr::vector<graph::key>, std::pmr::vector<colour>>
          &sorted_neighs,
      std::pmr::vector<graph::key> &clique, std::size_t &max_clique_size,
      std::size_t upper_bound, std::size_t &num_nodes,
      lru_cache<std::pmr::vector<graph::key>, std::pmr::vector<colour>> &cache,
      std::size_t &cache_hits, std::size_t depth = 1);

  std::pair<std::pmr::vector<graph::key>, std::pmr::vector<colour>> colour_sort(
      const std::pmr::vector<graph::key> &neighs,
      lru_cache<std::pmr::vector<graph::key>, std::pmr::vector<colour>> &cache,
      std::size_t &cache_hits) const {
    static thread_local size_t cache_hits_mod = 10'000;
    auto precomputed = cache.get(neighs);
    if (precomputed.has_value()) {
      cache_hits++;
      if (cache_hits % cache_hits_mod == 0) {
        cache_hits_mod += .5 * cache_hits_mod;
        logger::warn("current cache hits", cache_hits);
      }
      return precomputed.value();
    } else {
      auto [_, colours] = G.colour_sort(neighs);
      return cache.set(neighs, std::move(colours));
    }
  }

public:
  static inline std::size_t no_upper_bound = -1u;

  explicit multithreaded(const graph &g) : G(g), max_clique(memory::pool()) {
    logger::info("Number of available Threads", args::num_threads);
  }

  ~multithreaded() { logger::debug("~multithreaded()"); }

  std::pmr::vector<vertex>
  solve(flavour algo = flavour::exact,
        // the expected behaviour is (as far as I have tested) the function call
        // with be launched with the up-to-date values, even though the it seems
        // to be at compile time, it is dynamic initialisation
        std::size_t lower_bound = args::lower_bound,
        std::size_t upper_bound = args::upper_bound);

  void draw(const std::pmr::vector<vertex> &clique) { // E.draw(clique);
  }
};
} // namespace mc
#endif // MAXCLIQUE_HPP
