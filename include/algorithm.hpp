#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "graph.hpp"

namespace mc {
struct thread_id {
  thread_id(std::uint16_t id) : id(id) {}
  std::uint16_t id;
};

struct clique : thread_id {
  clique(thread_id self = default_thread, size_type size = 0)
      : thread_id(self), heuristic_size(size) {}

  graph::vertex_set nodes;
  size_type heuristic_size;
};

struct source_clique : clique {
  source_clique() = default;

  bool state_owner_unsafe(const thread_id &target) const {
    return id == target.id;
  }

  bool state_owner(const thread_id &target) {
    std::scoped_lock state_lck(state_mtx);
    return id == target.id;
  }

  size_type heuristic_unsafe() const { return heuristic_size; }

  size_type heuristic() {
    std::scoped_lock state_lck(state_mtx);
    return heuristic_size;
  }

  bool heuristic(clique &target, size_type new_size) {
    assert(new_size > target.heuristic_size);
    if (std::shared_lock state_shared_lck(state_mtx);
        new_size > heuristic_size) {
      state_shared_lck.unlock();
      {
        std::unique_lock state_unique_lck(state_mtx);
        const size_type old_size = heuristic_size;
        if ((heuristic_size = std::max(new_size, old_size)) > old_size) {
          id = target.id; // set state owner
        }
      }
      state_shared_lck.lock();
      target.heuristic_size = heuristic_size;
      return id == target.id;
    } else {
      return false;
    }
  }

  void flush(clique &target) {
    if (std::scoped_lock clique_lck(state_mtx);
        target.nodes.size() > nodes.size() && state_owner_unsafe(target)) {
      heuristic_size = target.nodes.size();
      nodes = std::move(target.nodes);
    }
  }

  std::shared_mutex state_mtx;
};

struct sync_branching {
  sync_branching(const graph &g, source_clique &source)
      : g(g), pruned(g.num_vertices()), source(source) {
    vertices.reserve(g.num_vertices());
    std::iota(vertices.begin(), vertices.end(), 0);
  }

  const graph &g;
  std::vector<bool> pruned;
  source_clique &source;
  graph::vertex_set vertices;
  std::atomic_bool abort_search = false, upper_bound_reached = false;

  std::mutex fetch_mtx;
  size_type upper_bound = -1u;

  graph::vertex_set fetch(graph::vertex v,
                          std::unique_lock<std::mutex> &fetch_lck) {
    auto all_neighs = g.neighbours(v);
    graph::vertex_set neighs;
    neighs.reserve(all_neighs.size());
    if (not fetch_lck.owns_lock()) {
      fetch_lck.lock();
    }
    std::copy_if(std::execution::par, all_neighs.begin(), all_neighs.end(),
                 std::back_inserter(neighs),
                 [this](graph::vertex u) { return not pruned[u]; });
    pruned[v] = true;
    fetch_lck.unlock();
    if (neighs.empty()) {
      abort_search = true;
    }

    return neighs;
  }
};

struct bound {
  graph::vertex_set &vertices;

  operator bool() { return not reached(); }

  bool reached() { return colours.empty(); }

  graph::vertex next() {
    assert(not vertices.empty());
    graph::vertex u = vertices.back();
    vertices.pop_back();
    return u;
  }
};

using colour = size_type;
struct colouring_bound {
  sync_branching &sync;
  const clique &local;
  std::vector<colour> colours;

  colouring_bound(graph::vertex_set &vertices, sync_branching &sync,
                  const clique &local)
      : vertices(vertices), sync(sync), local(local) {}

  operator bool() { return not reached(); }

  bool reached() { return colours.empty(); }

  graph::vertex next() {
    assert(not vertices.empty());
    graph::vertex u = vertices.back();
    vertices.pop_back();
    return u;
  }

  bool operator()() {
    std::vector<graph::vertex_set> colour_classes(vertices.size());
    size_type num_colours = 0;
    for (const auto &current : vertices) {
      colour col = 0;
      while (std::any_of(colour_classes[col].begin(), colour_classes[col].end(),
                         [&current, this](graph::vertex u) {
                           return sync.g.adjacent(u, current);
                         })) {
        col++;
      }
      colour_classes[col].emplace_back(current);
      num_colours = std::max(num_colours, col + 1);
    }

    if (num_colours <= local.heuristic_size) {
      return false;
    }

    colours.reserve(vertices.size());
    vertices.clear();
    for (colour col = 0; col < num_colours; ++col) {
      std::fill_n(std::back_inserter(colours), colour_classes[col].size(),
                  col + 1);
      std::copy(colour_classes[col].begin(), colour_classes[col].end(),
                std::back_inserter(vertices));
    }

    return true;
  }
};

struct bound_strategy {
  bound_strategy(std::unique_lock<std::mutex> &fetch_lck, graph::vertex v,
                 sync_branching &sync, clique &local)
      : sync(sync), local(local), v(v),
        heuristic_orig(update_local_heuristic()),
        vertices(sync.fetch(v, fetch_lck)) {}

  bound_strategy(bound_strategy &other)
      : sync(other.sync), local(other.local), v(other.next()),
        depth(other.depth + 1), heuristic_orig(update_local_heuristic()),
        vertices(sync.g.neighbourhood(v, other.vertices)) {}

  ~bound_strategy() {
    if (heuristic_orig < local.heuristic_size &&
        sync.source.state_owner_unsafe(local)) {
      local.nodes.emplace_back(v);
    }
  }

  size_type update_local_heuristic() {
    return local.heuristic_size = sync.source.heuristic();
  }

  bool update_heuristic(size_type new_size) {
    return new_size > local.heuristic_size
               ? sync.source.heuristic(local, new_size)
               : false;
  }

  bool prepare() {
    if (not vertices.empty() || depth == sync.upper_bound ||
        sync.upper_bound_reached) {
      if (update_heuristic(depth)) {
        if (depth == sync.upper_bound) {
          sync.upper_bound_reached = true;
        }
        local.nodes.clear();
        local.nodes.reserve(depth);
      }
      return false;
    } else {
      return initialise();
    }
  }

  virtual operator bool() { return not reached(); }

  virtual bool reached() { return vertices.empty(); }

  virtual graph::vertex next() {
    assert(not vertices.empty());
    graph::vertex u = vertices.back();
    vertices.pop_back();
    return u;
  }

protected:
  virtual bool initialise() = 0;

  sync_branching &sync;
  clique &local;

  const graph::vertex v;
  const size_type depth = 1;
  const size_type heuristic_orig;
  graph::vertex_set vertices;
};

struct bound_by_colouring : bound_strategy {
  bound_by_colouring(std::unique_lock<std::mutex> &fetch_lck, graph::vertex v,
                     sync_branching &sync, clique &local)
      : bound_strategy(fetch_lck, v, sync, local) {}

  bound_by_colouring(bound_by_colouring &other) : bound_strategy(other) {}

  bool reached() final {
    return bound_strategy::reached() || colours.back() <= local.heuristic_size;
  }

  graph::vertex next() final {
    graph::vertex v = bound_strategy::next();
    colours.pop_back();
    return v;
  }

protected:
  bool initialise() final {
    std::vector<graph::vertex_set> colour_classes(vertices.size());
    size_type num_colours = 0;
    for (const auto &current : vertices) {
      colour col = 0;
      while (std::any_of(colour_classes[col].begin(), colour_classes[col].end(),
                         [&current, this](graph::vertex u) {
                           return sync.g.adjacent(u, current);
                         })) {
        col++;
      }
      colour_classes[col].emplace_back(current);
      num_colours = std::max(num_colours, col + 1);
    }

    if (num_colours <= local.heuristic_size) {
      return false;
    }

    colours.reserve(vertices.size());
    vertices.clear();
    for (colour col = 0; col < num_colours; ++col) {
      std::fill_n(std::back_inserter(colours), colour_classes[col].size(),
                  col + 1);
      std::copy(colour_classes[col].begin(), colour_classes[col].end(),
                std::back_inserter(vertices));
    }

    return true;
  }

  std::vector<colour> colours;
};

// template <typename BoundStrategy>
using BoundStrategy = bound_by_colouring;
struct branch_clique : clique {
  static_assert(std::is_base_of_v<bound_strategy, BoundStrategy>);

  using bound_strategy_type = BoundStrategy;
  branch_clique(thread_id self, sync_branching &sync)
      : clique(self, sync.source.heuristic()), sync(sync),
        fetch_lck(sync.fetch_mtx, std::try_to_lock) {}

  ~branch_clique() { sync.source.flush(*this); }

  virtual void run(bound_strategy_type &bound) = 0;

  bool operator()(graph::vertex v) {
    bound_strategy_type bound(fetch_lck, v, sync, *this);
    if (bound.prepare()) {
      run(bound);
      return true;
    } else {
      return false;
    }
  }

  sync_branching &sync;
  std::unique_lock<std::mutex> fetch_lck;
};

struct algorithm_heuristic : branch_clique {
  void run(bound_strategy_type &bound) final {
    bound_strategy_type next(bound);
    if (next.prepare()) {
      run(next);
    }
  }
};

struct algorithm_exact : branch_clique {
  void run(bound_strategy_type &bound) final {
    while (bound) {
      bound_strategy_type next(bound);
      if (next.prepare()) {
        run(next);
      }
    }
  }
};

// template <typename Algorithm> struct algorithm_impl : source_clique {
//   static_assert(
//       std::is_base_of_v<typename Algorithm::bound_strategy, bound_strategy>);

//   algorithm_impl(const graph &g)
//       : source_clique(g), threads(num_threads) {}

//   graph::vertex_set operator()() {

//     while (not vertices.empty() && not predicate()) {
//       {
//         std::unique_lock thread_lck(thread_mtx);
//         pool.wait(thread_lck, [this] {
//           return std::any_of(available.begin(), available.end(),
//                              [](bool th) { return th; });
//         });
//       }
//       selector.prepare();
//       for (std::uint16_t th_id = 0; th_id < threads.size() && not
//       predicate();
//            ++th_id) {
//         if (available[th_id]) {
//           available[th_id] = false;
//           threads[th_id] = selector.next();
//         }
//       }
//     }

//     return std::move(nodes);
//   }

//   std::shared_mutex thread_mtx;
//   std::condition_variable_any pool;
//   std::vector<std::thread> threads;
//   std::vector<bool> available;
// };

// template <typename Algorithm> struct algorithm_hybrid {};

// struct algo_naive : algorithm_impl {
//   void setup(const graph &g) final;
//   void cleanup() final;
// };

// struct algo_degeneracy : algorithm_impl {
//   void setup(const graph &) final {}
//   void cleanup() final {}
// };

// struct algo_colours : algorithm_impl {
//   void setup(const graph &) final {}
//   void cleanup() final {}
// };

// using algorithm = std::variant<algo_naive, algo_degeneracy, algo_colours>;
} // namespace mc

#endif
