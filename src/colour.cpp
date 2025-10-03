#include "colour.hpp"

namespace mc {

colouring_workspace::colouring_workspace(const graph &G)
    : G(G), V(G.vertex_count()), colours(V, memory::pool()),
      degrees(V, memory::pool()), epoch_mark(V, memory::pool()),
      bucket_pos(V, memory::pool()), bucket_head(V, key_npos, memory::pool()),
      next(V, key_npos, memory::pool()), prev(V, key_npos, memory::pool()),
      out_keys(memory::pool()), out_colours(memory::pool()),
      colour_mask(memory::pool()), touched(memory::pool()),
      order(memory::pool()) {}

void colouring_workspace::order_smallest_last(std::pmr::vector<key> &R) {
  degree min_degree = R.size();

  for (key v : R) {
    epoch_mark[v] = epoch;
  }

  for (key v : R) {
    degree d = 0;
    for (auto [u, last] = G.neighbours(v); u != last; ++u) {
      d += (epoch_mark[*u] == epoch);
    }

    degrees[v] = d;
    bucket_pos[v] = d;

    next[v] = bucket_head[d];
    if (bucket_head[d] != key_npos) {
      prev[bucket_head[d]] = v;
    }
    bucket_head[d] = v;

    if (d < min_degree) {
      min_degree = d;
    }
  }

  std::size_t order_idx = 0;

  for (std::size_t k = 0; k < R.size(); ++k) {
    while (min_degree <= R.size() && bucket_head[min_degree] == key_npos) {
      ++min_degree;
    }

    key v = bucket_head[min_degree];
    bucket_head[min_degree] = next[v];
    if (next[v] != key_npos) {
      prev[next[v]] = key_npos;
    }

    R[order_idx++] = v;
    epoch_mark[v] = 0;

    for (auto [u, last] = G.neighbours(v); u != last; ++u) {
      if (epoch_mark[*u] != epoch) {
        continue;
      }

      degree old_d = degrees[*u];
      degree new_d = old_d - 1;

      if (prev[*u] != key_npos) {
        next[prev[*u]] = next[*u];
      } else {
        bucket_head[old_d] = next[*u];
      }

      if (next[*u] != key_npos) {
        prev[next[*u]] = prev[*u];
      }

      next[*u] = bucket_head[new_d];
      if (bucket_head[new_d] != key_npos) {
        prev[bucket_head[new_d]] = *u;
      }
      prev[*u] = key_npos;
      bucket_head[new_d] = *u;

      degrees[*u] = new_d;
      bucket_pos[*u] = new_d;
      if (new_d < min_degree) {
        min_degree = new_d;
      }
    }
  }

  std::reverse(R.begin(), R.end());
}

void colouring_workspace::greedy_colour(const std::pmr::vector<key> &R) {
  for (key v : R) {
    touched.clear();

    for (auto [first, last] = G.neighbours(v); first != last; ++first) {
      key u = *first;
      if (epoch_mark[u] == epoch) {
        colour cu = colours[u];
        if (cu > 0) {
          const std::size_t w = static_cast<std::size_t>(cu) >> 6;
          if (w >= colour_mask.size()) {
            colour_mask.resize(w + 1, 0);
          }
          uint64_t b = uint64_t(1) << (cu & 63);
          if ((colour_mask[w] & b) == 0) {
            colour_mask[w] |= b;
            touched.emplace_back(cu);
          }
        }
      }
    }

    colour c = first_free_colour();
    colours[v] = c;
    epoch_mark[v] = epoch;
    out_keys.emplace_back(v);
    out_colours.emplace_back(c);

    for (colour tc : touched) {
      colour_mask[tc >> 6] &= ~(uint64_t(1) << (tc & 63));
    }
  }
}

bool colouring_workspace::valid_colouring(const colour_sorted &cs) const {
  std::pmr::unordered_map<key, colour> colour_of(memory::pool());
  colour_of.reserve(cs.size());
  for (size_t i = 0; i < cs.size(); ++i) {
    colour_of[cs.key_at(i)] = cs.colour_at(i);
  }
  for (key u : cs.get_keys()) {
    auto [nb, nb_end] = G.neighbours(u);
    for (auto it = nb; it != nb_end; ++it) {
      key v = *it;
      if (u < v) {
        auto itc = colour_of.find(v);
        if (itc != colour_of.end() && itc->second == colour_of[u]) {
          return false;
        }
      }
    }
  }
  return true;
}
} // namespace mc
