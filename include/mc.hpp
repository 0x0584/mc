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

#include "branch.hpp"

namespace mc {
class clique_solver {
  const graph &G;
  branch::global_context ctx;

public:
  static inline const std::size_t no_upper_bound = -1u;

  explicit clique_solver(const graph &g) : G(g) {
    logger::info("Number of available Threads", args::num_threads);
  }

  std::pmr::vector<graph::vertex>
  solve(flavour algo = flavour::exact,
        std::size_t lower_bound = args::lower_bound,
        std::size_t upper_bound = args::upper_bound);

private:
  void solution(flavour algo, std::size_t upper_bound,
                const colour_sorted &root);

  inline bool explore_branch(key k, colour c) {
    if (ctx.bounded.load(std::memory_order_acquire)) [[unlikely]] {
      return false;
    } else if (c <= ctx.snapshot_size()) [[unlikely]] {
      logger::debug(G.to_vertex(k), "has no sufficient colours. abort search.");
      return false;
    } else {
      return true;
    }
  }

  static inline void set_intersection(const key __restrict *nf,
                                      const key __restrict *nl,
                                      const key __restrict *ef,
                                      const key __restrict *el,
                                      std::pmr::vector<key> &out) {
    while (nf != nl && ef != el) {
      if (*nf < *ef) {
        out.emplace_back(*nf++);
      } else if (*ef < *nf) {
        ++ef;
      } else {
        ++nf;
        ++ef;
      }
    }
    std::copy(nf, nl, std::back_inserter(out));
  }

  inline bool induce_neighbours(key k, colour_sorted root, std::size_t i,
                                std::pmr::vector<key> &neighs) {
    auto [nf, nl] = G.neighbours(k);
    neighs.reserve(G.degree(k));

    std::ptrdiff_t at = static_cast<std::ptrdiff_t>(i);
    std::pmr::vector<key> excluded(root.get_keys().begin() + at,
                                   root.get_keys().end(), memory::pool());

    std::sort(excluded.begin(), excluded.end());
    auto ef = excluded.data();
    set_intersection(nf, nl, ef, ef + excluded.size(), neighs);

    if (neighs.empty()) [[unlikely]] {
      logger::debug(G.to_vertex(k), "has no neighbours. abandon branching.");
      return false;
    } else {
      return true;
    }
  }
};

} // namespace mc
#endif // MAXCLIQUE_HPP
