#ifndef FAST_MAXCLIQUE_HPP
#define FAST_MAXCLIQUE_HPP

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <climits>
#include <cstring>
#include <execution>
#include <functional>
#include <iomanip>
#include <iostream>
#include <istream>
#include <map>
#include <memory_resource>
#include <mutex>
#include <numeric>
#include <set>
#include <shared_mutex>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#define COL_GREEN "\x1b[32m"
#define COL_MAGENTA "\x1b[35m"
#define COL_CYAN "\x1b[36m"
#define COL_RESET "\x1b[0m"
#define COL_BOLD "\x1b[1m"

//#define LOG_SLEEP
#ifndef LOG
#undef LOG_SLEEP
#endif

#ifdef PROFILER
#include "profiler.h"
#endif

#ifndef THREADS_PER_CORE
#define THREADS_PER_CORE 8
#endif

//#define NDEBUG

#define LITERAL(expr) #expr

#define mc_range(cont) (cont).begin(), (cont).end()
#define mc_const_range(cont) (cont).cbegin(), (cont).cend()

// #define LOG

using namespace std::chrono_literals;

struct log {
  static inline const auto short_hold = 1s, long_hold = 3s;

#ifdef LOG
  static inline std::mutex mtx;
#endif

  enum flags { none = 0x0, bold = 0b0001, ansi_colours = 0b0010 };

  static void setup_logger() {
    std::cout << std::fixed << std::setprecision(3) << std::left;
    std::cerr << std::fixed << std::setprecision(3) << std::left;
  }

  template <typename Chrono>
  static inline double duration(Chrono begin, Chrono end) {
    return std::chrono::duration<double>(end - begin).count();
  }

  template <typename Chrono>
  static inline std::string time_diff(Chrono begin, Chrono end,
                                      int flags = ansi_colours) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << std::left;
    if (flags & ansi_colours) {
      oss << COL_GREEN;
    }
    if (flags & bold) {
      oss << COL_BOLD;
    }
    oss << duration(begin, end) << "s";
    if (flags & ansi_colours) {
      oss << COL_RESET;
    }
    return oss.str();
  };

  static inline std::string progress(const std::size_t index,
                                     const std::size_t size) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << std::left
        << ((double(index + 1) * 100 / size)) << "%";
    return oss.str();
  }

  template <typename Chrono> static inline void pause(Chrono time) {
#ifdef LOG_SLEEP
    std::this_thread::sleep_for(time);
#else
    (void)time;
#endif
  }

  template <typename... Args> static inline void info(Args &&...args) {
#ifdef LOG
    std::scoped_lock print_lock(mtx);
#endif
    std::cerr << "> ";
    ((std::cerr << std::forward<Args>(args) << " "), ...);
    std::cerr << COL_RESET << std::endl;
  }

  template <typename... Args> static inline void print(Args &&...args) {
#ifdef LOG
    std::scoped_lock print_lock(mtx);
#endif
    ((std::cout << std::forward<Args>(args) << " "), ...);
    std::cout << COL_RESET << std::endl;
  }

  template <typename... Args>
  static inline void print_thread(std::uint32_t thread_id, Args &&...args) {
#ifdef LOG
    std::scoped_lock print_lock(mtx);
#endif
    std::cout << std::string((thread_id + 1), ' ') << std::to_string(thread_id)
              << " ";
    ((std::cout << std::forward<Args>(args) << " "), ...);
    std::cout << COL_RESET << std::endl;
  }
};

namespace thread {
const std::uint32_t threads_per_core = THREADS_PER_CORE;
const std::uint32_t num_threads =
    std::thread::hardware_concurrency() * threads_per_core;

void set_priority(std::thread &th, int policy, int priority);

template <typename Callable> struct scope_dtor {
  scope_dtor(const scope_dtor &) = delete;
  scope_dtor(scope_dtor &&) = delete;

  inline explicit scope_dtor(Callable &&fn)
      : callback(std::forward<Callable>(fn)) {}
  inline ~scope_dtor() { callback(); }

  scope_dtor &operator=(const scope_dtor &) = delete;
  scope_dtor &operator=(scope_dtor &&) = delete;

private:
  Callable callback;
};

}; // namespace thread

namespace mc {
using size_t = std::size_t;

template <typename T>
static inline void container_reserve_memory(T &container, std::size_t size) {
  container.max_load_factor(0.5);
  container.reserve(size);
}

class input {
  std::size_t num_v, num_e, size = -1u;
  bool expect_size;
  std::istream &in;

  static inline bool fetch(std::istringstream &iss, std::size_t &value) {
    if (iss.bad() || iss.eof()) {
      return false;
    } else {
      iss >> value;
      return true;
    }
  }

public:
  // TODO: handle header info
  inline input(std::istream &in) : in(in) {
    std::string line;
    std::getline(in, line);
    std::istringstream iss(line);
    if (not(fetch(iss, num_v) && fetch(iss, num_e))) {
      throw std::runtime_error("Could not fetch input header");
    }
    if ((expect_size = fetch(iss, size))) {
      log::info("Expected Max Clique of size", size);
    }
  }

  inline std::istream &operator*() { return in; }
  inline std::istream *operator->() { return &in; }

  inline std::size_t num_vertices() const { return num_v; }
  inline std::size_t num_edges() const { return num_e; }
  inline std::size_t max_clique_size() const { return size; }

  inline bool is_maximal_size(const std::size_t clique_size) const {
    return not expect_size || clique_size >= size;
  }
};

struct graph {
  friend struct enumerator;

  using vertex = unsigned;
  using neighbours_t = std::unordered_set<vertex>;
  using adjacency_t = std::unordered_map<vertex, std::unordered_set<vertex>>;

  graph(const graph &) = default;
  graph(graph &&) = default;

  inline graph(input &in, bool undirected = true)
      : in(in), undirected(undirected) {
    auto begin = std::chrono::high_resolution_clock::now();

    container_reserve_memory(adj_lst, in.num_vertices());
    for (std::size_t n_edges = in.num_edges();
         not in->eof() && in->good() && n_edges--;) {
      vertex u, v;
      *in >> u >> v; // TODO: handle comments
      add_edge(u, v);
    }

    auto end = std::chrono::high_resolution_clock::now();
    log::info((undirected ? "Undirected" : "Directed"), "Graph with",
              in.num_vertices(), "vertices and", edge_count,
              "edges was read in", log::time_diff(begin, end, log::bold));
    log::pause(log::short_hold);

    assert(adj_lst.size() == in.num_vertices());
  }

  inline void add_edge(const vertex u, const vertex v) {
    if (vertex_alloc(u); edge_add(u, v)) {
      ++edge_count;
      if (vertex_alloc(v); undirected) {
        edge_add(v, u);
      }
    }
  }

  inline const neighbours_t &neighbours(const vertex v) const {
    assert(adj_lst.count(v));
    return adj_lst.at(v);
  }

  inline const adjacency_t &adjacency_list() const { return adj_lst; }

  void print() const;

private:
  inline void vertex_alloc(const vertex u) {
    if (adj_lst[u].empty()) {
      container_reserve_memory(adj_lst[u], in.num_vertices());
    }
  }

  inline bool edge_add(const vertex u, const vertex v) {
    return adj_lst[u].emplace(v).second;
  }

  const input &in;
  const bool undirected;

  mc::size_t edge_count = 0;
  adjacency_t adj_lst;
};

// vertex wrapper that serves as a medium to access vertices as keys, in order
// to take advantage of std::vector (since when executing the algorithm, the
// set of vertices would remain constant) instead of how they are stored as
// std::unordered_map within the graph, it acts also as a handler of vertices
// for colouring and inducing vertex-neighbourhood
struct enumerator {
  using adjacency_t =
      std::vector<std::pair<graph::vertex, graph::neighbours_t>>;

  using key = unsigned;
  using colour = unsigned;

  const input &in;

  enumerator(graph &&g);

  inline graph::vertex vertex(const mc::size_t index) const {
    assert(index < count());
    return adj_lst_orig[index].first;
  }

  inline std::vector<key> neighbourhood(const key v,
                                        const std::vector<key> &neighs) const {
    assert(v < count());
    std::vector<key> new_neighs;
    new_neighs.reserve(neighs.size());
    for (const key u : neighs) {
      if (adj_mat[v][u]) { // neighbours of both vertices u and v
        new_neighs.emplace_back(u);
      }
    }
    return new_neighs;
  }

  inline std::vector<key>
  neighbours(const key v, std::function<bool(const key)> &&probe) const {
    assert(v < count());
    std::vector<key> neighs;
    neighs.reserve(adj_lst[v].size());
    for (const key u : adj_lst[v]) {
      if (probe(u)) {
        neighs.emplace_back(u);
      }
    }
    return neighs;
  }

  inline std::vector<graph::vertex>
  unfold_keys(const std::vector<key> &keys) const {
    std::vector<graph::vertex> vertices(keys.size());
    std::transform(std::execution::par_unseq, mc_const_range(keys),
                   vertices.begin(), [this](const key v) {
                     assert(v < count());
                     return vertex(v);
                   });

    std::set<graph::vertex> m(mc_const_range(vertices));
    std::ostringstream oss;
    oss << "Max Clique has " << m.size() << " vertices { ";
    for (const auto &e : m) {
      oss << e << " ";
    }
    oss << "}";
    log::info(oss.str());

    return vertices;
  }

  std::vector<colour> greedy_colour_sort(std::vector<key> &neighs) const;

  bool is_clique(const std::vector<key> &clique) const;

  inline std::size_t count() const { return in.num_vertices(); }

  void print() const;

private:
  adjacency_t adj_lst_orig;
  std::vector<std::vector<key>>
      adj_lst; // Adjacency List for fast neighbourhood deduction
  std::vector<std::vector<bool>>
      adj_mat; // Adjacency Matrix for fast edge probing
};

class fast_maxclique {
public:
  enum struct algorithm_type { exact, heuristic, hybrid };

  static inline const mc::size_t maximum_bound = -1u;

private:
  const enumerator &enumtor;

  std::atomic_bool upper_bound_reached = false;

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
  mc::size_t overall_size = 0;
  //
  // after branch termination, if the current thread had found the largest one
  // so far amongst all the running threads (even if they are still running
  std::vector<enumerator::key> max_clique;
  //
  // hence, we can set the max clique few times and avoid unnecessary
  //  assignments of cliques from several threads, at least in most cases
  std::uint32_t holder_thread_id = -1u;

  void solution(algorithm_type algo, mc::size_t upper_bound);

  inline bool enlarge_clique_size(const std::uint32_t thread_id,
                                  mc::size_t &max_clique_size,
                                  const mc::size_t depth);

  void branch_exact(const std::uint32_t thread_id, const enumerator::key v,
                    std::vector<enumerator::key> &neighs,
                    std::vector<enumerator::colour> &colours,
                    std::vector<enumerator::key> &clique,
                    mc::size_t &max_clique_size, const mc::size_t upper_bound,
                    std::size_t &num_nodes, const mc::size_t depth = 1);

  void branch_heuristic(const std::uint32_t thread_id, const enumerator::key v,
                        std::vector<enumerator::key> &neighs,
                        std::vector<enumerator::key> &clique,
                        mc::size_t &max_clique_size,
                        const mc::size_t upper_bound, std::size_t &num_nodes,
                        const mc::size_t depth = 1);

public:
  static inline const mc::size_t no_upper_bound = -1u;

  fast_maxclique(const enumerator &enumtor) : enumtor(enumtor) {
    log::pause(log::long_hold);
  }

  friend std::ostream &operator<<(std::ostream &os, algorithm_type algo_type) {
    if (algo_type == algorithm_type::exact) {
      return os << "Exact Algorithm";
    } else if (algo_type == algorithm_type::heuristic) {
      return os << "Heuristic Algorithm";
    } else {
      return os << "Hybrid Algorithm";
    }
  }

  std::vector<graph::vertex>
  find_maxclique(algorithm_type algo = algorithm_type::exact,
                 mc::size_t lower_bound = 0u,
                 mc::size_t upper_bound = no_upper_bound);
};
} // namespace mc

namespace thread {
struct pool {};
} // namespace thread
#endif
