#ifndef BRANCH_HPP
#define BRANCH_HPP

#include <atomic>
#include <shared_mutex>
#include <vector>

#include "colour.hpp"

namespace mc {
namespace branch {
struct global_context {
  std::shared_mutex mtx;
  std::size_t global_size = 0;
  std::pmr::vector<key> global_clique{memory::pool()};
  key owner_key = key();
  std::atomic_bool bounded = false;

  inline bool update_max_clique(key root, std::size_t depth) {
    if (std::shared_lock r(mtx); depth <= global_size) {
      return false;
    }
    if (std::scoped_lock w(mtx); depth > global_size) {
      global_size = depth;
      owner_key = root;
      logger::debug("Enlarged clique size to", global_size);
      return true;
    }
    return false;
  }

  inline std::size_t snapshot_size() {
    std::shared_lock r(mtx);
    return global_size;
  }

  inline void reset() {
    std::scoped_lock w(mtx);
    global_clique.clear();
    global_size = 0;
    owner_key = key();
    bounded.store(false, std::memory_order_release);
  }
};

#ifndef NDEBUG
struct stats {
  key root;
  std::size_t visited_nodes = 0;
  std::size_t expanded_nodes = 0;
  std::size_t pruned_by_colour = 0;
  std::size_t pruned_by_chromatic = 0;
  std::size_t max_depth_reached = 0;
  std::size_t cumulative_depth = 0;

  const char *tag;

  explicit stats(key r, const char *policy) : root(r), tag(policy) {}

  inline void on_visit(std::size_t depth) noexcept {
    thread_local auto last_log = std::chrono::steady_clock::now();
    ++visited_nodes;
    cumulative_depth += depth;
    if (depth > max_depth_reached) {
      max_depth_reached = depth;
    }

    auto now = std::chrono::steady_clock::now();
    if (now - last_log >= std::chrono::seconds(5)) {
      log_summary();
      last_log = now;
    }
  }

  inline void on_expand() noexcept { ++expanded_nodes; }
  inline void on_prune_colour() noexcept { ++pruned_by_colour; }
  inline void on_prune_chromatic() noexcept { ++pruned_by_chromatic; }

  inline void log_summary() const {
    const std::size_t total_pruned = pruned_by_colour + pruned_by_chromatic;

    const double d_visited_nodes = static_cast<double>(visited_nodes);
    const double prune_ratio =
        visited_nodes ? total_pruned / d_visited_nodes : 0.0;
    const double branch_factor =
        visited_nodes ? expanded_nodes / d_visited_nodes : 0.0;
    const double avg_depth =
        visited_nodes ? cumulative_depth / d_visited_nodes : 0.0;

    std::ostringstream oss;
    oss << tag << " root=" << root << " visited=" << visited_nodes
        << " expanded=" << expanded_nodes
        << " pruned_colour=" << pruned_by_colour
        << " pruned_chromatic=" << pruned_by_chromatic
        << " pruned_ratio=" << std::fixed << std::setprecision(3)
        << (prune_ratio * 100.0) << "%"
        << " branch_factor=" << std::setprecision(3) << branch_factor
        << " max_depth=" << max_depth_reached
        << " avg_depth=" << std::setprecision(3) << avg_depth;

    logger::debug(oss.str());
  }
};
#endif

namespace policy {
struct exact {
  static constexpr const char *tag() noexcept { return "Exact"; }

  inline colour_sorted make_child(colouring_engine &engine,
                                  std::pmr::vector<key> &&keys) const {
    return engine.colour_sort(std::move(keys));
  }

  inline void descend(std::pmr::vector<colour_sorted> &parents,
                      std::pmr::vector<key> &path, key u, colour_sorted &&child,
                      colour_sorted &) const {
    path.emplace_back(u);
    parents.emplace_back(std::move(child));
  }

  inline void on_leaf_backtrack(std::pmr::vector<colour_sorted> &,
                                std::pmr::vector<key> &) const {
    // no-op: keep peeling siblings at this depth
  }

  inline void on_prune_chromatic(std::pmr::vector<colour_sorted> &,
                                 std::pmr::vector<key> &) const {
    // no-op: keep peeling siblings
  }
};

struct heuristic {
  static constexpr const char *tag() noexcept { return "Heuristic"; }

  inline colour_sorted make_child(colouring_engine &engine,
                                  std::pmr::vector<key> &&keys) const {
    return engine.colour_sort(std::move(keys));
  }

  inline void descend(std::pmr::vector<colour_sorted> &,
                      std::pmr::vector<key> &path, key u, colour_sorted &&child,
                      colour_sorted &parent) const {
    path.emplace_back(u);
    parent = std::move(child); // overwrite at same depth
  }

  inline void on_leaf_backtrack(std::pmr::vector<colour_sorted> &parents,
                                std::pmr::vector<key> &path) const {
    parents.pop_back();
    path.pop_back();
  }

  inline void on_prune_chromatic(std::pmr::vector<colour_sorted> &parents,
                                 std::pmr::vector<key> &path) const {
    parents.pop_back();
    path.pop_back();
  }
};
} // namespace policy

template <typename Policy> class runner {
  const graph &G;
  key root;
  branch::global_context &ctx;
  colouring_engine &engine;
  const Policy &policy;

  std::pmr::vector<colour_sorted> parents{memory::pool()};
  std::pmr::vector<key> path{memory::pool()};

public:
  std::pmr::vector<key> clique{memory::pool()};

#ifndef NDEBUG
  branch::stats stats;
#endif

  runner(const graph &g, key r, branch::global_context &ctx,
         colouring_engine &engine, const Policy &policy = Policy{})
      : G(g), root(r), ctx(ctx), engine(engine), policy(policy)
#ifndef NDEBUG
        ,
        stats(r, Policy::tag())
#endif
  {
  }

  ~runner() {
#ifndef NDEBUG
    stats.log_summary();
#endif
  }

  bool find_clique(colour_sorted R, std::size_t upper_bound) {
    bool found_new_clique = run(std::move(R), upper_bound);
    if (root == ctx.owner_key && found_new_clique) {
      std::scoped_lock lock(ctx.mtx);
      if (clique.size() == ctx.global_size) {
        ctx.global_clique = std::move(clique);
        if constexpr (logger::current_level == logger::log_level::debug) {
          std::ostringstream oss;
          oss << "Found clique for " << G.to_vertex(root) << " of "
              << ctx.global_clique.size() << " vertices { ";
          for (key u : ctx.global_clique)
            oss << G.to_vertex(u) << " ";
          oss << "}";
          logger::debug(oss.str());
        }
      }
    }
    return found_new_clique;
  }

private:
  bool run(colour_sorted &&R, std::size_t upper_bound) {
    prepare(std::move(R));

    while (!parents.empty() && !ctx.bounded.load(std::memory_order_acquire))
        [[likely]] {
      auto &parent = parents.back();
      if (parent.empty()) [[unlikely]] {
        backtrack_parent();
        continue;
      }

      auto [u, c] = parent.peel();
      std::size_t depth = path.size();

#ifndef NDEBUG
      stats.on_visit(depth);
#endif

      if (colour_bound(depth, c)) {
        continue;
      }

      auto child_keys = G.neighbourhood(u, parent.get_keys());

      depth++;
      if (depth == upper_bound) [[unlikely]] {
        ctx.bounded.store(true, std::memory_order_release);
      }

      if (child_keys.empty() || ctx.bounded.load(std::memory_order_acquire)) {
        if (commit_leaf(depth, u)) {
          return true;
        } else {
          continue;
        }
      }

      auto child = policy.make_child(engine, std::move(child_keys));
      if (chromatic_hit(depth, child)) {
        continue;
      }

      descend(u, std::move(child), parent);
    }

    return false;
  }

  inline void backtrack_parent() {
    parents.pop_back();
    path.pop_back();
  }

  inline void prepare(colour_sorted &&R) {
    parents.clear();
    path.clear();
    clique.clear();
    parents.emplace_back(std::move(R));
    path.emplace_back(root);
  }

  inline bool colour_bound(std::size_t depth, colour c) {
    if (depth + c > ctx.snapshot_size()) [[likely]] {
      return false;
    } else {
#ifndef NDEBUG
      stats.on_prune_colour();
#endif
      backtrack_parent();
      return true;
    }
  }

  inline bool commit_leaf(std::size_t leaf_depth, key u) {
    if (!ctx.update_max_clique(root, leaf_depth)) {
      policy.on_leaf_backtrack(parents, path);
      return false;
    }

    clique.clear();
    clique.reserve(leaf_depth);
    clique.insert(clique.end(), path.begin(), path.end());
    clique.emplace_back(u);

#ifndef NDEBUG
    stats.log_summary();
#endif

    return true;
  }

  inline bool chromatic_hit(std::size_t leaf_depth,
                            const colour_sorted &child) {
    Assert(child.empty(), "this should never trigger!");
    const std::size_t clique_potential_size =
        leaf_depth + child.chromatic_num();
    if (clique_potential_size > ctx.snapshot_size()) [[unlikely]] {
      return false;
    } else {
      policy.on_prune_chromatic(parents, path);
#ifndef NDEBUG
      stats.on_prune_chromatic();
#endif
      return true;
    }
  }

  inline void descend(key u, colour_sorted &&child, colour_sorted &parent) {
    policy.descend(parents, path, u, std::move(child), parent);
#ifndef NDEBUG
    stats.on_expand();
#endif
  }
};
} // namespace branch
} // namespace mc

#endif
