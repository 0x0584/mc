#ifndef COLOUR_HPP
#define COLOUR_HPP

#include "graph.hpp"

#include <numeric>
#include <stdexcept>

namespace mc {

class colour_sorted {
  std::pmr::vector<key> keys;
  std::pmr::vector<colour> colours;

public:
  colour_sorted() : keys(memory::pool()), colours(memory::pool()) {}

  inline void reserve(std::size_t sz) {
    keys.reserve(sz);
    colours.reserve(sz);
  }

  inline const std::pmr::vector<key> &get_keys() const { return keys; }

  inline const std::pmr::vector<colour> &get_colours() const { return colours; }

  inline bool empty() const noexcept { return colours.empty(); }

  inline std::size_t size() const noexcept { return colours.size(); }

  inline void clear() {
    keys.clear();
    colours.clear();
  }

  inline void emplace_back(key u, colour c) {
    keys.emplace_back(u);
    colours.emplace_back(c);
  }

  inline colour chromatic_num() const {
    if (empty())
      throw std::runtime_error("chromatic_num() on empty colour_sorted");
    return colours.back();
  }

  inline std::pair<key, colour> peel() {
    if (empty())
      throw std::runtime_error("peel() on empty colour_sorted");
    std::pair<key, colour> kc{keys.back(), colours.back()};
    keys.pop_back();
    colours.pop_back();
    return kc;
  }

  inline std::pair<key, colour> at(std::size_t i) const {
    if (empty())
      throw std::runtime_error("at() on empty colour_sorted");
    return {keys[i], colours[i]};
  }

  inline key key_at(std::size_t i) const {
    if (empty())
      throw std::runtime_error("key_at() on empty colour_sorted");
    return keys[i];
  }

  inline key colour_at(std::size_t i) const {
    if (empty())
      throw std::runtime_error("colour_at() on empty colour_sorted");
    return colours[i];
  }
};

// FIXME: remove workspace
class colouring_workspace {
  const graph &G;
  const std::size_t V;
  timestamp epoch = 1;

  std::pmr::vector<colour> colours;
  std::pmr::vector<degree> degrees;
  std::pmr::vector<timestamp> epoch_mark;

  std::pmr::vector<key> bucket_pos;
  std::pmr::vector<key> bucket_head;
  std::pmr::vector<key> next;
  std::pmr::vector<key> prev;

  std::pmr::vector<key> out_keys;
  std::pmr::vector<key> out_colours;
  std::pmr::vector<uint64_t> colour_mask;
  std::pmr::vector<colour> touched;
  std::pmr::vector<key> order;

public:
  colouring_workspace() = delete;
  colouring_workspace(const colouring_workspace &) = delete;
  colouring_workspace &operator=(const colouring_workspace &) = delete;
  colouring_workspace(colouring_workspace &&) = default;
  colouring_workspace &operator=(colouring_workspace &&) = delete;

  colouring_workspace(const graph &G);

  inline void start_epoch(std::size_t sz) noexcept {
    ++epoch;

    colour_mask.clear();
    touched.clear();

    out_keys.clear();
    out_keys.reserve(sz);

    out_colours.clear();
    out_colours.reserve(sz);

    order.resize(sz);
    std::iota(order.begin(), order.end(), 0u);
  }

  inline void colour_sort(std::pmr::vector<key> &&R, colour_sorted &out,
                          bool order_keys) {
    if (R.empty()) [[unlikely]] {
      return;
    }

    start_epoch(R.size());

    if (order_keys) {
      order_smallest_last(R);
    }

    greedy_colour(R);

    stable_sort_by_colour(out);

    Assert(valid_colouring(out), "invalid colouring");
  }

private:
  void order_smallest_last(std::pmr::vector<key> &R);

  void greedy_colour(const std::pmr::vector<key> &R);

  inline colour first_free_colour() const {
    for (std::size_t w = 0; w < colour_mask.size(); ++w) {
      uint64_t inv = ~colour_mask[w];
      if (inv) {
        for (unsigned b = 0; b < 64; ++b) {
          if (inv & (static_cast<uint64_t>(1) << b)) {
            colour candidate = b + (w << 6);
            if (candidate > 0) {
              return candidate;
            }
          }
        }
      }
    }
    return 1u + (colour_mask.size() << 6);
  }

  inline void stable_sort_by_colour(colour_sorted &out) {
    std::stable_sort(order.begin(), order.end(), [&](key u, key v) {
      return out_colours[u] < out_colours[v];
    });
    for (key u : order) {
      out.emplace_back(out_keys[u], out_colours[u]);
    }
  }

  bool valid_colouring(const colour_sorted &cs) const;
};

class colouring_engine {
  colouring_workspace ws;

public:
  colouring_engine() = delete;
  colouring_engine(const colouring_engine &) = delete;
  colouring_engine &operator=(const colouring_engine &) = delete;
  colouring_engine(colouring_engine &&) = default;
  colouring_engine &operator=(colouring_engine &&) = delete;

  colouring_engine(const graph &G) : ws(G) {}

  static colour_sorted colour_sort(const graph &G) {
    std::pmr::vector<key> R(G.vertex_count(), memory::pool());
    std::iota(R.begin(), R.end(), 0u);
    colour_sorted out;
    out.reserve(R.size());
    // FIXME: optimise the recomputation of graph degrees
    colouring_workspace ws(G);
    ws.colour_sort(std::move(R), out, true);
    return out;
  }

  inline colour_sorted colour_sort_order_policy(std::pmr::vector<key> &&R,
                                                bool order_keys) {
    if (R.empty()) [[unlikely]] {
      return {};
    }

    colour_sorted out;
    out.reserve(R.size());
    ws.colour_sort(std::move(R), out, order_keys);

    return out;
  }

  inline colour_sorted colour_sort(std::pmr::vector<key> &&R) {
    return colour_sort_order_policy(std::move(R), true);
  }

  inline colour_sorted colour_sort_no_order(std::pmr::vector<key> &&R) {
    return colour_sort_order_policy(std::move(R), false);
  }
};

} // namespace mc
#endif
