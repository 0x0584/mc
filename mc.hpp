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

  static inline std::string progress(std::size_t index, std::size_t size) {
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

  template <typename... Args> static inline void printv(Args &&...args) {
    std::scoped_lock print_lock(mtx);
    ((std::cout << std::forward<Args>(args)), ...); // verbose
    std::cout << COL_RESET << std::endl;
  }

  template <typename... Args> static inline void error(Args &&...args) {
    std::scoped_lock print_lock(mtx);
    ((std::cerr << std::forward<Args>(args) << " "), ...);
    std::cerr << COL_RESET << std::endl;
    std::exit(EXIT_FAILURE);
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
    std::string filename;
    for (int ch; (ch = getopt(argc, argv, "r:i:s:u:l:deyh")) != -1;) {
      switch (ch) {
      case 'r': // XXX: handle --run N
        num_turns = std::max(1l, std::atol(optarg));
        break;
      case 'i': // XXX: handle --in input.g
        stdin = false;
        file = std::ifstream(filename = optarg);
        break;
      case 'd': // XXX: handle --edge-directed
        undirected = false;
        break;
      case 's': // XXX: handle --size N
        expect_size = true;
        size = static_cast<std::size_t>(std::atol(optarg));
        log::info("Expecting a Max Clique of size", size);
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
      case 'h':
      case ':':
      case '?':
      default:
        std::cerr
            << "Find the Max Clique of a graph using branch-and-bound\n\n"
            << "  -r N run for N times\n"
            << "  -i FILE take input from FILE instead of STDIN\n"
            << "  -u N expect a at most a clique of size N\n"
            << "  -l N expect at least a clique of size N\n"
            << "  -e run the algorithm as EXACT (default HEURISTIC)\n"
            << "  -y run the algorithm as HYBRID (HEURISTIC + EXACT)\n"
            << "  -d use DIRECTED edges instead of the default UNDIRECTED\n"
            << "\n";
        exit(EXIT_FAILURE);
      }
    }

    log::info("Running", args::exec_mode);

    if (stdin) {
      log::info("Reading from STDIN");
    } else {
      log::info("Reading from", filename);
    }
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
    if (not fetcher(num_v) || not fetcher(num_e)) {
      throw std::runtime_error("Could not parse input header, use -?");
    }
    if (not args::expect_size && (args::expect_size = fetcher(args::size))) {
      log::info("Expecting a Max Clique of size", args::size);
    }

    log::info("Source Graph is", (args::undirected ? "Undirected" : "Directed"),
              "with", num_v, "vertices and", num_e, "edges");
  }

  inline std::istream &operator*() { return args::stream(); }
  inline std::istream *operator->() { return &args::stream(); }

  inline bool is_maximal_size(std::size_t clique_size) const {
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
  // this is a primitive type, ins case of a change in the implementation
  // since some parts of the code should be updated to avoid overhead of
  // copying the  objects rather than either referencing them or moving them
  using vertex = unsigned;

  using neighbours_set = std::unordered_set<vertex>;
  using adjacency_map = std::unordered_map<vertex, neighbours_set>;

  friend struct enumerator;
  friend struct graph_builder;

  graph(const graph &) = default;
  graph(graph &&) = default;
  graph() = default;

  graph &operator=(graph &) = delete;
  graph &operator=(graph &&G) {
    this->A = std::move(G.A);
    this->edge_count = G.edge_count;
    this->undirected = G.undirected;
    return *this;
  }

  inline const neighbours_set &neighbours(vertex v) const {
    assert(A.count(v));
    return A.at(v);
  }

  inline const adjacency_map &adjacency() const { return A; }

  inline bool directed() { return not undirected; }

  void print() const {
    std::ostringstream oss;
    oss << (undirected ? "Undirected" : "Directed")
        << " Graph (vertices=" << A.size() << ", edges=" << edge_count << ")\n";
    for (const auto &[v, neighs] : A) {
      oss << v << " { ";
      for (vertex u : neighs) {
        oss << u << " ";
      }
      oss << "}\n";
    }
    oss << "\n";
    log::print(oss.str());
  }

private:
  bool undirected = false;
  mc::size_t edge_count = 0;
  adjacency_map A;
};

struct feed {
  static inline const std::int64_t CHUNK_SIZE = 16384, LINE_SIZE = 32;

  static inline const char deli = '\n', sep = ' ';

  using buffer = std::array<char, CHUNK_SIZE>;
  using line = std::array<char, LINE_SIZE>;
  using chunk = std::pair<buffer, std::size_t>;

  feed(feed &&feed) = delete;
  feed(const feed &feed) = delete;

  inline feed() : tail_remaining(remaining.begin()) {
    std::memset(remaining.data(), 0x0, remaining.size());
  }

  ~feed() {
    if (tail_remaining != remaining.begin()) {
      // buffer::difference_type size =
      //     std::distance(remaining.begin(), tail_remaining);
      // std::string s{remaining.data(), std::size_t(size)};
      // log::printv(" ==>>> `", s, "` with size=", size);
      log::error("INVALID file: no NL at the end of the file");
    }
  }

  inline operator bool() { return reading(); }
  inline std::size_t estimate_chunks() const {
    return 1 + in.num_e / CHUNK_SIZE;
  }
  inline std::size_t num_vertices() const { return in.num_v; }
  inline std::size_t num_edges() const { return in.num_e; }

  chunk read_chunk() {
    buffer buff;
    const buffer::difference_type size_remaining =
        std::distance(remaining.begin(), tail_remaining);
    std::move(remaining.begin(), tail_remaining, buff.data());

    auto read = read_next(CHUNK_SIZE - size_remaining);
    auto begin_read = read.first.begin();
    const auto size_read = read.second;
    std::move(begin_read, begin_read + size_read, buff.data() + size_remaining);

    const buffer::iterator tail_buff =
        buff.begin() + (size_remaining + size_read);
    buffer::iterator delimiter = tail_buff;
    while (delimiter != buff.begin() && *--delimiter != feed::deli)
      ;
    if (delimiter == buff.begin()) {
      log::error("FAILURE: chunk_size=", CHUNK_SIZE,
                 " exceeded! recompile with a bigger size");
    }

    chunk chunk{std::move(buff), std::distance(buff.begin(), delimiter++)};
    std::move(delimiter, tail_buff, remaining.data());
    tail_remaining = remaining.begin() + std::distance(delimiter, tail_buff);
    // log::printv(reads++, " remaining=`",
    //             std::string{remaining.begin(),
    //                         std::size_t(std::distance(delimiter,
    //                         tail_buff))},
    //             "` EOF=", in->eof());
    return chunk;
  }

private:
  inline bool reading() {
    if (in->bad()) {
      log::error("UNEXPECTED READ FAILURE");
    }
    return not in->eof() && not in->fail();
  }

  inline std::pair<buffer, const std::streamsize>
  read_next(std::streamsize size_read) {
    buffer read;
    if (in->read(read.data(), size_read); in->bad()) {
      log::error("FAILURE: cannot read from stream");
    }
    return {read, in->gcount()};
  }

  input in;
  buffer remaining{};
  buffer::iterator tail_remaining;
};

struct graph_builder {
  graph_builder(graph_builder &&) = delete;
  graph_builder(const graph_builder &) = delete;

  explicit graph_builder(bool undirected = args::undirected) {
    G.undirected = undirected;
    container_reserve_memory(G.A, feed.num_vertices());
    Q.reserve(feed.estimate_chunks());
  }

  graph_builder &operator=(const graph_builder &) = delete;
  graph_builder &operator=(graph_builder &&) = delete;

  bool read_single_vertex(graph::vertex &w, std::string::iterator &it,
                          std::string::iterator end) {
    auto skip = [end](std::string::iterator &it) {
      while (it != end && (*it == feed::sep || *it == feed::deli)) {
        ++it;
      }
    };
    auto read_vertex = [end](std::string::iterator &it) {
      int limit = 0;
      while (limit++ < 12 && it != end && *it != feed::sep &&
             *it != feed::deli) {
        ++it;
      }
    };

    if (skip(it); it != end) {
      std::string::iterator i = it;
      read_vertex(it);
      w = std::atol(std::string{i, it}.c_str());
      return true;
    } else {
      return false;
    }
  }

  graph build() {
    // TODO: improve reading vertices by finding optimal way to at the same time
    // knmow how many vertices are there and add edges to the graph

    auto begin = std::chrono::high_resolution_clock::now();

    std::deque<feed::chunk> q;
    std::unordered_map<graph::vertex, std::size_t> adj_count;
    container_reserve_memory(adj_count, feed.num_vertices());
    std::mutex adj_count_mtx;

    do {
      feed::chunk chunk = feed.read_chunk();
      std::string buffer{chunk.first.data(), chunk.second};
      q.push_back(std::move(chunk));
      Q.emplace_back(std::async(
          std::launch::deferred,
          [this, &adj_count, &adj_count_mtx](std::string buffer) {
            // it it useless to save vertices since we do not know how much
            // space we need to allocate anyway, that is why we count them first
            std::string::iterator it = buffer.begin();
            std::unordered_map<graph::vertex, std::size_t> adj_count_local;
            int edge_count = 2;
            for (graph::vertex u, v; edge_count == 2;) {
              edge_count = read_single_vertex(u, it, buffer.end()) +
                           read_single_vertex(v, it, buffer.end());
              if (edge_count != 2) {
                break;
              }
              ++adj_count_local[u];
              if (G.undirected) {
                ++adj_count_local[v];
              }
            }
            if (edge_count != 0) {
              log::error("could not read an edge", edge_count);
            }
            std::unique_lock<std::mutex> lock(adj_count_mtx);
            // to be benchmarked which style is bettter
            adj_count.insert(adj_count_local.begin(), adj_count_local.end());
          },
          std::move(buffer)));
    } while (feed);

    std::for_each(Q.begin(), Q.end(), [](std::future<void> &f) { f.get(); });
    Q.clear();

    for (auto [v, size] : adj_count) {
      container_reserve_memory(G.A[v], size);
    }

    std::atomic_size_t edge_count{0};
    while (not q.empty()) {
      feed::chunk &chunk = q.back();
      std::string buffer{std::move(chunk.first.data()), chunk.second};
      q.pop_back();
      Q.emplace_back(std::async(
          std::launch::deferred,
          [this, &edge_count](std::string buffer) {
            std::string::iterator it = buffer.begin();
            for (graph::vertex u, v; read_single_vertex(u, it, buffer.end()) &&
                                     read_single_vertex(v, it, buffer.end());) {
              std::unique_lock lock(graph_mtx);
              if (G.A[u].emplace(v).second) {
                ++edge_count;
                if (G.undirected) {
                  G.A[v].emplace(u);
                }
              }
            }
          },
          std::move(buffer)));
    }
    std::for_each(Q.begin(), Q.end(), [](std::future<void> &f) { f.get(); });

    G.edge_count = edge_count;
    auto end = std::chrono::high_resolution_clock::now();

    log::info("Graph with", feed.num_vertices(), "vertices and", G.edge_count,
              "edges was read in", log::time_diff(begin, end, log::bold));

    assert(G.A.size() == feed.num_vertices());
    assert(G.edge_count == feed.num_edges());

    LONG_DELAY();

    return std::move(G);
  }

private:
  std::mutex graph_mtx;
  std::vector<std::future<void>> Q;
  mc::feed feed;
  graph G;
};

// vertex wrapper that serves as a medium to access vertices as keys, in order
// to take advantage of std::vector (since when executing the algorithm, the
// set of vertices would remain constant) instead of how they are stored as
// std::unordered_map within the graph, it acts also as a handler of vertices
// for colouring and inducing vertex-neighbourhood
struct enumerator {
  // this is a premitive type too, same as graph::vertex so changes in the
  // implementation are required in order to avoid overhead of copying
  // instead of using references or moving the object
  using key = unsigned;
  using colour = unsigned;
  using adjacency_vector =
      std::vector<std::pair<graph::vertex, graph::neighbours_set>>;

public:
  enumerator();

  inline graph::vertex key_to_vertex(mc::size_t index) const {
    assert(index < vertex_count());
    return E[index].first;
  }

  inline std::vector<key> neighbourhood(key v,
                                        const std::vector<key> &neighs) const {
    assert(v < vertex_count());
    std::vector<key> new_neighs;
    new_neighs.reserve(neighs.size());
    std::copy_if(neighs.begin(), neighs.end(), std::back_inserter(new_neighs),
                 // neighbours of both vertices u and v
                 [this, v](key u) { return B[v][u]; });
    return new_neighs;
  }

  inline std::vector<key> neighbours(key v,
                                     std::unordered_set<key> &pruned) const {
    assert(v < vertex_count());
    std::vector<key> neighs;
    neighs.reserve(A.at(v).size());
    std::copy_if(A.at(v).begin(), A.at(v).end(), std::back_inserter(neighs),
                 [&pruned](key u) { return not pruned.count(u); });
    return neighs;
  }

  inline std::vector<graph::vertex>
  unfold_keys(const std::vector<key> &keys) const {
    std::vector<graph::vertex> vertices(keys.size());
    std::transform(std::execution::par_unseq, keys.begin(), keys.end(),
                   vertices.begin(), [this](key v) {
                     assert(v < vertex_count());
                     return key_to_vertex(v);
                   });
    return vertices;
  }

  void print() const {
    std::ostringstream oss;
    oss << "Enumerated Vertices\n";
    for (const auto &[v, neighs] : E) {
      oss << v << " { ";
      for (graph::vertex u : neighs) {
        oss << u << " ";
      }
      oss << "}\n";
    }
    oss << "\n";
    log::print(oss.str());
  }

  std::vector<colour> greedy_colour_sort(std::vector<key> &neighs) const;

  bool is_clique(const std::vector<key> &clique) const;

  inline std::size_t vertex_count() const { return E.size(); }

private:
  adjacency_vector E; // Enumertaed vertices

  template <typename T> using vector_2d = std::vector<std::vector<T>>;
  vector_2d<key> A;  // Adjacency List for fast neighbourhood deduction
  vector_2d<bool> B; // Adjacency Matrix for fast edge probing
};

class multithreaded {
public:
  static inline const mc::size_t maximum_bound = -1u;

  enumerator g;

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

  bool enlarge_clique_size(std::uint32_t thread_id, mc::size_t &max_clique_size,
                           mc::size_t depth);

  void branch_exact(std::uint32_t thread_id, enumerator::key v,
                    std::vector<enumerator::key> &neighs,
                    std::vector<enumerator::colour> &colours,
                    std::vector<enumerator::key> &clique,
                    mc::size_t &max_clique_size, mc::size_t upper_bound,
                    std::size_t &num_nodes, mc::size_t depth = 1);

  void branch_heuristic(std::uint32_t thread_id, enumerator::key v,
                        std::vector<enumerator::key> &neighs,
                        std::vector<enumerator::key> &clique,
                        mc::size_t &max_clique_size, mc::size_t upper_bound,
                        std::size_t &num_nodes, mc::size_t depth = 1);

public:
  static inline mc::size_t no_upper_bound = -1u;

  multithreaded() {
    log::info("Number of available Threads", thread::num_threads);
    LONG_DELAY();
  }

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
