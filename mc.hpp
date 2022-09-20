#ifndef MAXCLIQUE_HPP
#define MAXCLIQUE_HPP

#define LITERAL(expr) #expr
#define EMPTY_MACRO                                                            \
  do {                                                                         \
  } while (false)

#ifndef LOG
#define DELAY() EMPTY_MACRO
#define LONG_DELAY() EMPTY_MACRO
#else
#ifndef PAUSE_LONG_DELAY
#define PAUSE_LONG_DELAY 3s
#endif // PAUSE_DELAY_LONG
#ifndef PAUSE_DELAY
#define PAUSE_DELAY 2s
#endif // PAUSE_DELAY
#define LONG_DELAY() std::this_thread::sleep_for(PAUSE_LONG_DELAY)
#define DELAY() std::this_thread::sleep_for(PAUSE_DELAY)
#endif // LOG

#include <atomic>
#include <chrono>
#include <execution>
#include <functional>
#include <memory_resource>
#include <mutex>
#include <shared_mutex>
#include <thread>
#include <utility>

using namespace std::chrono_literals;

#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#define mc_range(cont) (cont).begin(), (cont).end()
#define mc_const_range(cont) (cont).cbegin(), (cont).cend()

//#define NDEBUG

#include <cassert>
#include <climits>
#include <cstring>
#include <getopt.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <sstream>
#include <string>

#define COL_GREEN "\x1b[32m"
#define COL_MAGENTA "\x1b[35m"
#define COL_CYAN "\x1b[36m"
#define COL_RESET "\x1b[0m"
#define COL_BOLD "\x1b[1m"

struct log {
  enum flags { none = 0x0, bold = 0b0001, ansi_colours = 0b0010 };

  static inline void setup_logger() {
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
  }

  static inline std::string progress(const std::size_t index,
                                     const std::size_t size) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << std::left
        << ((double(index + 1) * 100 / size)) << "%";
    return oss.str();
  }

  template <typename... Args> static inline void info(Args &&...args) {
    std::scoped_lock print_lock(mtx);
    std::cerr << "> ";
    ((std::cerr << std::forward<Args>(args) << " "), ...);
    std::cerr << COL_RESET << std::endl;
  }

  template <typename... Args> static inline void print(Args &&...args) {
    std::scoped_lock print_lock(mtx);
    ((std::cout << std::forward<Args>(args) << " "), ...);
    std::cout << COL_RESET << std::endl;
  }

  template <typename... Args>
  static inline void print_thread(std::uint32_t thread_id, Args &&...args) {
    std::scoped_lock print_lock(mtx);
    std::cout << std::string((thread_id + 1), ' ') << std::to_string(thread_id)
              << " ";
    ((std::cout << std::forward<Args>(args) << " "), ...);
    std::cout << COL_RESET << std::endl;
  }

private:
  static inline std::mutex mtx;
};

#ifndef THREADS_PER_CORE
#define THREADS_PER_CORE 8
#endif

namespace thread {
const std::uint32_t threads_per_core = THREADS_PER_CORE;
const std::uint32_t num_threads =
    std::thread::hardware_concurrency() * threads_per_core;

// since the thread can return on several conditions, it is
// practical to use a scope destructor to ensure that threads are
// marked as available not matter the branch
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

static inline void set_priority(std::thread &th, int policy, int priority) {
  sched_param sch_params;
  sch_params.sched_priority = priority;
  if (pthread_setschedparam(th.native_handle(), policy, &sch_params)) {
    log::info("Failed to set Thread scheduling:", std::strerror(errno));
    exit(-1);
  }
}
}; // namespace thread

namespace mc {
using size_t = std::size_t;

template <typename T>
static inline void container_reserve_memory(T &container, std::size_t size) {
  container.max_load_factor(0.5);
  container.reserve(size);
}

enum struct flavour { exact, heuristic, hybrid };

std::ostream &operator<<(std::ostream &os, flavour algo_type) {
  if (algo_type == flavour::exact) {
    return os << "Exact Algorithm";
  } else if (algo_type == flavour::heuristic) {
    return os << "Heuristic Algorithm";
  } else {
    return os << "Hybrid Algorithm";
  }
}

struct args {
  static void parse(int argc, char *argv[]) {
    for (int ch; (ch = getopt(argc, argv, "r:i:s:u:l:dey")) != -1;) {
      switch (ch) {
      case 'r': // XXX: handle --run N
        num_turns = std::max(1l, std::atol(optarg));
        break;
      case 'i': // XXX: handle --in input.g
        stdin = false;
        file = std::ifstream(optarg);
        break;
      case 'd': // XXX: handle --edge-directed
        undirected = false;
        break;
      case 's': // XXX: handle --size N
        expect_size = true;
        size = static_cast<std::size_t>(std::atol(optarg));
        log::info("Expected Max Clique of size", size);
        break;
      case 'u': // XXX: handle --upper-bound N
        upper_bound = static_cast<std::size_t>(std::atol(optarg));
        log::info("Expected Upper Bound for Max Clique of size", upper_bound);
        break;
      case 'l': // XXX: handle --lower-bound N
        lower_bound = static_cast<std::size_t>(std::atol(optarg));
        log::info("Expected Lower Bound for Max Clique of size", lower_bound);
        break;
      case 'e': // XXX: handle --exact
        if (exec_mode == flavour::heuristic) {
          // just in case both hybrid and excat were specified, run as hybrid
          exec_mode = flavour::exact;
        }
        break;
      case 'y': // XXX: handle --hybrid
        exec_mode = flavour::hybrid;
        break;
      case ':':
      case '?':
      default:
        std::cerr
            << "Find the Max Clique of a graph using branch-and-bound\n\n"
            << "  -r N run for N times\n"
            << "  -i FILE take input from FILE instead of STDIN\n"
            << "  -u N expect a at most a clique of size N\n"
            << "  -l N expect at least a clique of size N\n"
            << "  -e run the algorithm as EXACT branching (default HEURISTIC)\n"
            << "  -y run the algorithm as HYBRID (HEURISTIC + EXACT)"
            << "  -d use DIRECTED edges instead of the default UNDIRECTED\n";
        exit(EXIT_FAILURE);
      }
    }

    log::info("Running", exec_mode);
  }

  static inline std::istream &stream() { return stdin ? std::cin : file; }

  static inline std::ifstream file;
  static inline long num_turns = 100;
  static inline bool expect_size, undirected = true, stdin = true;
  static inline std::size_t size = -1u, upper_bound = -1u, lower_bound = 1;

  static inline flavour exec_mode = flavour::heuristic;
};

struct input {
  // TODO: handle header info
  //   V        |                 E                                  | S
  //   vertices | [un]directed edges (default undirected if omitted) | size
  inline input() {
    std::string line;
    std::getline(args::stream(), line);
    fetch_ftor fetcher(line);
    if (not(fetcher(num_v) && fetcher(num_e))) {
      throw std::runtime_error("Could not fetcher input header, use -?");
    }
    if (not args::expect_size && (args::expect_size = fetcher(args::size))) {
      log::info("Expected Max Clique of size", args::size);
    }

    log::info("Source Graph is", (args::undirected ? "Undirected" : "Directed"),
              "with", num_v, "vertices and", num_e, "edges");
  }

  inline std::istream &operator*() { return args::stream(); }
  inline std::istream *operator->() { return &args::stream(); }

  inline bool is_maximal_size(const std::size_t clique_size) const {
    return not args::expect_size || clique_size >= args::size;
  }

  struct fetch_ftor {
    fetch_ftor(std::string &line) : iss(line) {}

    template <typename T> bool operator()(T &value) {
      if (iss.bad() || iss.eof()) {
        return false;
      } else {
        iss >> value;
        return true;
      }
    }

  private:
    std::istringstream iss;
  };

  std::size_t num_v, num_e;
};

struct graph {
  using vertex = unsigned;

  using neighbours_set = std::unordered_set<vertex>;
  using adjacency_map = std::unordered_map<vertex, neighbours_set>;

  friend struct enumerator;

  graph(const graph &) = default;
  graph(graph &&) = default;

  inline graph(input &in) : graph(in, args::undirected) {}
  inline graph(input &in, bool undirected) : in(in), undirected(undirected) {
    auto begin = std::chrono::high_resolution_clock::now();

    container_reserve_memory(adj_lst, in.num_v);
    for (std::size_t n_edges = in.num_e;
         not in->eof() && in->good() && n_edges--;) {
      std::string line;
      // TODO: if extended file source is to be supported, optimised read must
      // be multithreaded by reading chunks and merging results of the running
      // threads based on the number of lines, i.e. number of edges defined
      if (std::getline(args::stream(), line);
          // XXX: handle comments when line starts with %
          not line.empty() && line.front() != '%') {
        vertex u, v;
        input::fetch_ftor fetcher = line;
        if (fetcher(u) && fetcher(v)) {
          add_edge(u, v);
        }
      }
      // log::info("Input left", log::progress(n_edges, in.num_e));
    }

    auto end = std::chrono::high_resolution_clock::now();
    log::info("Graph with", in.num_v, "vertices and", edge_count,
              "edges was read in", log::time_diff(begin, end, log::bold));
    assert(adj_lst.size() == in.num_v);

    LONG_DELAY();
  }

  inline void add_edge(const vertex u, const vertex v) {
    if (vertex_alloc(u); edge_add(u, v)) {
      ++edge_count;
      if (vertex_alloc(v); undirected) {
        edge_add(v, u);
      }
    }
  }

  inline const neighbours_set &neighbours(const vertex v) const {
    assert(adj_lst.count(v));
    return adj_lst.at(v);
  }

  inline const adjacency_map &adjacency_list() const { return adj_lst; }

  void print() const {
    std::ostringstream oss;
    oss << (undirected ? "Undirected" : "Directed")
        << " Graph (vertices=" << adj_lst.size() << ", edges=" << edge_count
        << ")\n";
    for (const auto &[v, neighs] : adj_lst) {
      oss << v << " { ";
      for (const auto &u : neighs) {
        oss << u << " ";
      }
      oss << "}\n";
    }
    oss << "\n";
    log::print(oss.str());
  }

private:
  inline void vertex_alloc(const vertex u) {
    if (adj_lst[u].empty()) {
      container_reserve_memory(adj_lst[u], in.num_v);
    }
  }

  inline bool edge_add(const vertex u, const vertex v) {
    return adj_lst[u].emplace(v).second;
  }

  const input &in;
  const bool undirected;

  mc::size_t edge_count = 0;
  adjacency_map adj_lst;
};

// vertex wrapper that serves as a medium to access vertices as keys, in order
// to take advantage of std::vector (since when executing the algorithm, the
// set of vertices would remain constant) instead of how they are stored as
// std::unordered_map within the graph, it acts also as a handler of vertices
// for colouring and inducing vertex-neighbourhood
struct enumerator {
  using key = unsigned;
  using colour = unsigned;
  using adjacency_vector =
      std::vector<std::pair<graph::vertex, graph::neighbours_set>>;

  input in;
  graph g;

public:
  enumerator();

  inline graph::vertex key_to_vertex(const mc::size_t index) const {
    assert(index < vertex_count());
    return adj_lst_orig[index].first;
  }

  inline std::vector<key> neighbourhood(const key v,
                                        const std::vector<key> &neighs) const {
    assert(v < vertex_count());
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
    assert(v < vertex_count());
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
                     assert(v < vertex_count());
                     return key_to_vertex(v);
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

  void print() const {
    std::ostringstream oss;
    oss << "Enumerated Vertices\n";
    for (const auto &[v, neighs] : adj_lst_orig) {
      oss << v << " { ";
      for (const auto &u : neighs) {
        oss << u << " ";
      }
      oss << "}\n";
    }
    oss << "\n";
    log::print(oss.str());
  }

  std::vector<colour> greedy_colour_sort(std::vector<key> &neighs) const;

  bool is_clique(const std::vector<key> &clique) const;

  inline std::size_t vertex_count() const { return in.num_v; }

private:
  adjacency_vector adj_lst_orig;

  template <typename T> using vector_2d = std::vector<std::vector<T>>;
  vector_2d<key> adj_lst;  // Adjacency List for fast neighbourhood deduction
  vector_2d<bool> adj_mat; // Adjacency Matrix for fast edge probing
};

class multithreaded {
public:
  static inline const mc::size_t maximum_bound = -1u;

private:
  enumerator g;

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

  // terminate the algorithm early if the depth matches the bound
  std::atomic_bool upper_bound_reached = false;

  void solution(flavour algo, mc::size_t upper_bound);

  bool enlarge_clique_size(const std::uint32_t thread_id,
                           mc::size_t &max_clique_size, const mc::size_t depth);

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
  static inline mc::size_t no_upper_bound = -1u;

  multithreaded() { LONG_DELAY(); }

  std::vector<graph::vertex>
  solve(flavour algo = flavour::exact,
        // the expected behaviour is (as far as I have tested) the function call
        // with be launched with the up-to-date values, even though the it seems
        // to be at compile time, it is dynamic initialisation
        mc::size_t lower_bound = args::lower_bound,
        mc::size_t upper_bound = args::upper_bound);
};
} // namespace mc
#endif // MAXCLIQUE_HPP
