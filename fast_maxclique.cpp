#include <cassert>
#include <pthread.h>
#include <cstring>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <execution>
#include <map>
#include <memory_resource>
#include <mutex>
#include <set>
#include <shared_mutex>
#include <sstream>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#define COL_GREEN    "\x1b[32m"
#define COL_MAGENTA  "\x1b[35m"
#define COL_CYAN     "\x1b[36m"
#define COL_RESET    "\x1b[0m"
#define COL_BOLD     "\x1b[1m"

#ifdef PROFILER
# include "profiler.h"
#endif

#ifndef THREADS_PER_CORE
# define THREADS_PER_CORE 8
#endif

//#define NDEBUG

using namespace std::chrono_literals;

namespace thread {
  void set_priority(std::thread &th, int policy, int priority) {
    sched_param sch_params;
    sch_params.sched_priority = priority;
    if (pthread_setschedparam(th.native_handle(), policy, &sch_params)) {
      std::cerr << "Failed to set Thread scheduling : " << std::strerror(errno) << std::endl;
      exit(-1);
    }
  }
}

#if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
# define INLINE __attribute__ ((always_inline)) inline
#elif defined(_MSC_VER) && _MSC_VER >= 1400
# define INLINE __forceinline inline
#else
# define INLINE inline
#endif

#define LITERAL(expr) #expr

template <typename Callable>
struct end_of_scope_executor {
  end_of_scope_executor(const end_of_scope_executor &) = delete;
  end_of_scope_executor(end_of_scope_executor &&) = delete;

    INLINE explicit end_of_scope_executor(Callable &&fn)
        : callback(std::forward<Callable>(fn)) {}
  INLINE ~end_of_scope_executor() { callback(); }

  end_of_scope_executor &operator=(const end_of_scope_executor &) = delete;
  end_of_scope_executor &operator=(end_of_scope_executor &&) = delete;

private:
  Callable callback;
};

INLINE double percentage(const std::size_t index, const std::size_t size) {
  return double(index + 1) * 100 / size;
}

namespace graph {
  using vertex = unsigned;
  using key = unsigned;
  using size_t = unsigned;
  using colour = unsigned;

  const std::uint32_t threads_per_core = THREADS_PER_CORE;
  const std::uint32_t num_threads = std::thread::hardware_concurrency() * threads_per_core;

  struct fast_maxclique {
    enum struct algorithms { exact, heuristic, hybrid };

    friend std::ostream &operator<<(std::ostream &os, algorithms algo_type) {
      if (algo_type == algorithms::exact) {
        return os << "Exact Algorithm";
      } else if (algo_type == algorithms::heuristic) {
        return os << "Heuristic Algorithm";
      } else {
        return os << "Hybrid Algorithm";
      }
    }

#define neighbours_at(index)                    adj_lst[index].second

    graph::size_t expected_max_clique_size, num_v, num_e;

    fast_maxclique() {
#ifdef PROFILER
      Profiler::beginSession("Preparing Graph", "load.json");
#endif

      auto begin = std::chrono::high_resolution_clock::now();

#ifdef PROFILER
      Profiler::InstrumentationTimer timer_reading("Reading from STDIN");
#endif
      std::cin >> num_v >> num_e >> expected_max_clique_size;
      std::unordered_map<graph::vertex, std::set<graph::vertex>> tmp_adj_lst;
      tmp_adj_lst.max_load_factor(0.5);
      tmp_adj_lst.reserve(num_v);
      for (size_t n_edges = num_e; n_edges--;) {
        graph::vertex u, v;
        std::cin >> u >> v;
        tmp_adj_lst[u].emplace(v);
        tmp_adj_lst[v].emplace(u);
      }
#ifdef PROFILER
      timer_reading.stop();
#endif

      assert(tmp_adj_lst.size() == num_v);

#ifdef PROFILER
      Profiler::InstrumentationTimer timer_sorting("Ordering vertices");
#endif

      std::unordered_map<graph::vertex, graph::size_t> neighbours_degree;
      neighbours_degree.max_load_factor(0.5);
      neighbours_degree.reserve(num_v);
      for (const auto &[v, neighbours] : tmp_adj_lst) {
        for (const auto &u : neighbours) {
          neighbours_degree[v] += tmp_adj_lst[u].size();
        }
      }
      std::vector<std::pair<graph::vertex, std::set<graph::vertex>>> adj_lst_orig;
      adj_lst_orig.reserve(num_v);
      for (auto &v : tmp_adj_lst) {
        adj_lst_orig.emplace_back(std::make_pair(v.first, std::move(v.second)));
      }
      std::sort(std::execution::par_unseq, adj_lst_orig.begin(), adj_lst_orig.end(),
                [&](const std::pair<graph::vertex, std::set<graph::vertex>> &v,
                    const std::pair<graph::vertex, std::set<graph::vertex>> &u) {
                  return (u.second.size() < v.second.size() ||
                          (u.second.size() == v.second.size() &&
                           neighbours_degree.at(u.first) < neighbours_degree.at(v.first)));
                });
#ifdef PROFILER
      timer_sorting.stop();
      Profiler::InstrumentationTimer timer_adj_lst("Creating Adjacency List");
#endif

      std::unordered_map<graph::vertex, graph::key> reversed_mapping;
      reversed_mapping.max_load_factor(0.5);
      reversed_mapping.reserve(num_v);

      adj_lst.reserve(num_v);
      adj_mtx.reserve(num_v);
      for (const auto &v : adj_lst_orig) {
        reversed_mapping.emplace(v.first, adj_lst.size());

        adj_mtx.emplace_back(std::vector<bool>(num_v));
        adj_lst.emplace_back(std::make_pair(v.first, std::vector<graph::key>{}));
        adj_lst.back().second.reserve(v.second.size());
      }

      for (const auto &[v, neighbours_orig] : adj_lst_orig) {
        const graph::key &v_key = reversed_mapping.at(v);
        std::vector<bool> &v_adj_mtx = adj_mtx[v_key];
        std::vector<graph::key> &neighbours = neighbours_at(v_key);
        for (const auto &u : neighbours_orig) {
          const graph::key &u_key = reversed_mapping.at(u);
          v_adj_mtx[u_key] = true;
          neighbours.emplace_back(u_key);
        }
        std::sort(neighbours.begin(), neighbours.end());
      }

#ifdef PROFILER
      timer_adj_lst.stop();
#endif

      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Graph (vertices=" << num_v << " edges=" << num_e
                << ") with expected maxlique=" << expected_max_clique_size
                << " was read in " << std::chrono::duration<double>(end-begin).count() << "s\n\n";
#ifdef PROFILER
      Profiler::endSession();
#endif

#if 0//def LOG
      for (const auto &v : adj_lst) {
        std::cout << "Vertex: " << v.first << " { ";
        for (const auto &u : v.second) {
          std::cout << adj_lst[u].first << " ";
        }
        std::cout << "}\n";
      }
      std::cout << "\n";

      std::map<graph::size_t, std::vector<graph::vertex>> degrees;
      for (const auto &v : adj_lst) {
        degrees[v.second.size()].emplace_back(v.first);
      }

      for (const auto &deg : degrees) {
        std::cerr << "Degree " << deg.first << ": { ";
        for (const auto &v : deg.second) {
          std::cerr << v << " ";
        }
        std::cerr << "}\n";
      }
      std::cerr << "\n";
#endif

#if 0//ndef LOG
      std::map<graph::size_t, std::vector<graph::vertex>> degrees;
      for (const auto &v : adj_lst) {
        degrees[v.second.size()].emplace_back(v.first);
      }
      for (const auto &deg : degrees) {
        std::cerr << "Degree " << deg.first << ": " << deg.second.size() << "\n";
      }
      std::cerr << "\n";
#endif
#ifdef LOG
      std::this_thread::sleep_for(3s);
#endif
    }

    static const graph::size_t no_upper_bound = -1u;
    std::vector<graph::vertex>
    find_maxclique(algorithms algo_type = algorithms::exact,
                   graph::size_t lower_bound = 0u, graph::size_t upper_bound = no_upper_bound) {
#ifdef PROFILER
      Profiler::beginSession("Max Clique", "solution.json");
#endif
      assert(lower_bound <= upper_bound);
      if (upper_bound > 1) {
        overall_max_clique_size = lower_bound > 0 ? lower_bound - 1 : 0;
        solution(algo_type, upper_bound);
      } else if (upper_bound == 1) {
        max_clique.emplace_back(0);
      }
#ifdef PROFILER
      Profiler::endSession();
#endif

      std::vector<graph::vertex> clique;
      clique.reserve(max_clique.size());
      for (const auto &v : max_clique) {
        clique.emplace_back(adj_lst[v].first);
      }

      std::set<graph::vertex> m(clique.begin(), clique.end());
      std::cout << "SIZE = " << m.size() << " { ";
      for (const auto &e : m) {
        std::cout << e << " ";
      }
      std::cout << "}\n\n";

      if (not is_clique(max_clique)) {
        throw std::runtime_error(" NOT A CLIQUE!!\n");
      }

      if (algo_type != algorithms::heuristic) {
        if (max_clique.size() < lower_bound) {
          throw std::runtime_error(" MAX CLIQUE LOWER BOUND HAVE NOT BEEN MET\n");
        } else if (upper_bound == no_upper_bound && expected_max_clique_size != max_clique.size()) {
          throw std::runtime_error(" NOT A MAX CLIQUE!!\n");
        }
      }
      max_clique.clear();

      overall_max_clique_size = 0;
      upper_bound_reached = false;
      owner_thread_id = -1u;

      return clique;
    }

  private:

    std::vector<graph::key> max_clique;
    std::vector<std::pair<graph::vertex, std::vector<graph::key>>> adj_lst; // Adjacency List
    std::vector<std::vector<bool>> adj_mtx; // Adjacency Matrix

#ifdef LOG
    mutable std::mutex print_mtx;
    void print(const std::string &s) const {
      std::scoped_lock print_lock(print_mtx);
      std::cerr << "> " << s << COL_RESET << std::endl;
    }

    void print(const std::string &s, std::uint32_t thread_id) const {
      std::scoped_lock print_lock(print_mtx);
      std::cerr << std::string((thread_id + 1), ' ') << std::to_string(thread_id)
                << " " << s << std::endl;
    }
#endif

    std::vector<graph::size_t> degrees() const {
      std::vector<graph::size_t> D;
      D.reserve(num_v);
      for (const auto &[_, neighbours] : adj_lst) {
        D.emplace_back(neighbours.size());
      }
      return D;
    }

    INLINE std::vector<graph::colour>
    sort_vertices_by_colour(std::vector<graph::key> &neighbours) const {
      if (neighbours.empty()) {
        return {};
      }

#ifdef PROFILER
      Profiler::InstrumentationTimer timer_thread(std::string("Colouring Vertices"));
#endif

      std::vector<std::vector<graph::key>> colour_classes;
      colour_classes.resize(neighbours.size());
      for (auto &colour_class : colour_classes) {
        colour_class.reserve(neighbours.size());
      }

      graph::size_t num_colours = 0;
      for (const auto &v : neighbours) {
        graph::colour colour = 0;
        while (std::any_of(colour_classes[colour].begin(), colour_classes[colour].end(),
                           [&v, this](const graph::key &u) { return adj_mtx[v][u]; })) {
          colour++;
        }
        colour_classes[colour].emplace_back(v);
        num_colours = std::max(num_colours, colour + 1);
      }

      std::vector<graph::colour> colours;
      colours.reserve(neighbours.size());

      neighbours.clear();

      for (graph::colour colour = 0; colour < num_colours; ++colour) {
        std::fill_n(std::back_inserter(colours), colour_classes[colour].size(), colour + 1);
        std::copy_n(colour_classes[colour].begin(), colour_classes[colour].size(), std::back_inserter(neighbours));
      }

      return colours;
    }

    std::shared_mutex max_clique_mtx;
    graph::size_t overall_max_clique_size = 0;
    std::atomic_bool upper_bound_reached = false;
    std::uint32_t owner_thread_id = -1u;

    template <typename Chrono>
    std::string time_diff_str(Chrono begin, Chrono end, bool ansi_colours = true) {
      std::ostringstream oss;
      oss << std::fixed << std::setprecision(3) << std::left;
      if (ansi_colours) {
        oss << COL_GREEN << COL_BOLD;
      }
      oss << std::chrono::duration<double>(end - begin).count() << "s";
      if (ansi_colours) {
        oss << COL_RESET;
      }
      return oss.str();
    }

    void solution(algorithms algo_type, graph::size_t upper_bound) {
      auto begin = std::chrono::high_resolution_clock::now();

      if (algo_type == algorithms::hybrid) {
        solution(algorithms::heuristic, upper_bound);
        if (upper_bound_reached) {
          auto end = std::chrono::high_resolution_clock::now();
          std::cout << algorithms::hybrid << " finished using " << algorithms::heuristic << " only! ("
                    << std::chrono::duration<double>(end - begin).count() << "s)\n";
          return;
        }
      }

      std::vector<graph::size_t> D = degrees();

      ///////////////////////

      std::shared_mutex thread_mtx;
      std::condition_variable_any thread_pool;
      std::vector<std::thread> threads(num_threads);
      std::vector<bool> available(num_threads, true);

      std::mutex read_mtx;
      std::condition_variable reading_is_complete;
      bool reading = false;

      std::atomic_size_t total_branches = 0;

      bool abort_search = false;
      owner_thread_id = -1u;

      /////////////////////////

#ifndef LOG
      graph::size_t old_max_clique_size = overall_max_clique_size;
#endif
#ifdef PROFILER
      Profiler::InstrumentationTimer timer_loop("Looping over vertices");
#endif
      std::vector<graph::key> vertices;
      vertices.reserve(num_v);
      for (graph::key v = 0; v < num_v; ++v) {
        vertices.emplace_back(v);
      }

      auto percent_begin = std::chrono::high_resolution_clock::now();
      while (not vertices.empty() && not abort_search && not upper_bound_reached) {
        auto percent_end = std::chrono::high_resolution_clock::now();

        const std::size_t percent = num_v - vertices.size();
#ifdef LOG
        {
          std::scoped_lock print_lock(print_mtx);
          std::cerr << "> waiting... " << percentage(percent, num_v) << "% in "
                    << time_diff_str(begin, percent_end) << "\n";
        }
#else
        if (std::chrono::duration<double>(percent_end - percent_begin).count() >= 1.) {
          std::cerr << percentage(percent, num_v) << "% in "
                    << time_diff_str(begin, percent_end, false) << "\n";
        }
#endif
        percent_begin = percent_end;

        {
          std::unique_lock thread_lock(thread_mtx);
          thread_pool.wait(thread_lock, [&] {
            return std::any_of(available.begin(), available.end(), [](bool th) { return th; });
          });
        }

        if (abort_search || upper_bound_reached) {
          break;
        }

        std::vector<graph::colour> colours = sort_vertices_by_colour(vertices);

        graph::size_t current_max_clique_size;
        {
          std::shared_lock clique_lock(max_clique_mtx);
          current_max_clique_size = overall_max_clique_size;
        }

#ifdef LOG
        print("new thread(s) available! max_clique=" + std::to_string(current_max_clique_size));
#else
        if (old_max_clique_size != current_max_clique_size) {
          std::cerr << " max_clique=" << current_max_clique_size << "\n";
          old_max_clique_size = current_max_clique_size;
        }
#endif

        for (std::uint32_t thread_id = 0; thread_id < threads.size() && not vertices.empty()
               && not abort_search && not upper_bound_reached; ++thread_id) {

          if (available[thread_id]) {
            available[thread_id] = false;
            if (threads[thread_id].joinable()) {
              threads[thread_id].join();
            }

            const graph::key v = vertices.back();
            if (colours.back() <= current_max_clique_size) {
              abort_search = true;
#ifdef LOG
              print(std::to_string(adj_lst[v].first) + " has insufficiant colours. abort search!");
#else
              std::cerr << "abort search! max_clique=" << current_max_clique_size << "\n";
#endif
              break;
            }


            vertices.pop_back();
            colours.pop_back();

            reading = true;
            threads[thread_id] = std::thread([&, this](const std::uint32_t thread_id, const graph::key v,
                                                       graph::size_t max_clique_size) {
              end_of_scope_executor set_thread_available([&available, &thread_id,
                                                          &thread_mtx, &thread_pool] {
                {
                  std::shared_lock thread_lock(thread_mtx);
                  available[thread_id] = true;
                }
                thread_pool.notify_one();
              });

#ifdef PROFILER
              Profiler::InstrumentationTimer timer_thread(std::string("Branching on ") + std::to_string(v));
              Profiler::InstrumentationTimer timer_vertex_reading(std::string("Reading neighbour vertices of ")
                                                                  + std::to_string(v));
#endif

              std::vector<graph::key> neighbours;
              neighbours.reserve(D[v]);

              for (const auto &u : neighbours_at(v)) {
                if (D[u] != 0) {
                  neighbours.emplace_back(u);
                }
              }
              D[v] = 0;

#ifdef PROFILER
              timer_vertex_reading.stop();
#endif

              {
                std::scoped_lock reading_lock(read_mtx);
                reading = false;
              }
              reading_is_complete.notify_one();

              if (neighbours.empty()) {
#ifdef LOG
                print(std::to_string(adj_lst[v].first) + " has no neighbours. abort search!");
#endif
                abort_search = true;
                return;
              }

              std::vector<graph::colour> colours = sort_vertices_by_colour(neighbours);
              if (colours.back() < max_clique_size) {
#ifdef LOG
                print(std::to_string(adj_lst[v].first) + " has insufficiant colours. abort search!");
#endif
                abort_search = true;

                return;
              }

#ifdef LOG
              print(std::string("branching on ") + std::to_string(adj_lst[v].first)
                    + " with " + std::to_string(neighbours.size()) + " neighbours has "
                    + std::to_string(colours.back()) + " colours",
                    thread_id);
              auto begin = std::chrono::high_resolution_clock::now();
#endif

              std::vector<graph::key> clique;
              std::size_t num_nodes = 0;
              if (algo_type == algorithms::heuristic) {
                branch_heuristic(thread_id, v, neighbours, clique, max_clique_size, upper_bound, num_nodes);
              } else {
                branch_exact(thread_id, v, neighbours, colours, clique, max_clique_size, upper_bound, num_nodes);
              }

              if (std::scoped_lock clique_lock(max_clique_mtx);
                  clique.size() > max_clique.size() && thread_id == owner_thread_id) {
#ifdef LOG
                std::ostringstream oss;
                oss << "found clique for " << adj_lst[v].first << " of size=" << clique.size() << " { ";
                for (const auto &u : clique) {
                  oss << adj_lst[u].first << " ";
                }
                oss << "}";
                print(oss.str(), thread_id);
#endif
                max_clique = std::move(clique);
                overall_max_clique_size = max_clique.size();
              }

              total_branches += num_nodes;

#ifdef LOG
              auto end = std::chrono::high_resolution_clock::now();
              {
                std::scoped_lock print_lock(print_mtx);
                std::cerr << std::string((thread_id + 1), ' ') <<  thread_id << " done with " << adj_lst[v].first
                          << " after " << num_nodes << " branches took "
                          << std::chrono::duration<double>(end - begin).count() << "s\n";
              }
#endif
            }, thread_id, v, current_max_clique_size);

            std::unique_lock reading_lock(read_mtx);
            reading_is_complete.wait(reading_lock, [&reading]() {
              return not reading;
            });

          }
        }
      }

#ifdef LOG
      print("stopped making threads");
#endif
      for (auto &th : threads) {
        if (th.joinable()) {
          th.join();
        }
      }
      auto end = std::chrono::high_resolution_clock::now();

#ifdef LOG
      std::cerr << "> all threads joined\n\n";
#endif
      std::cout << algo_type << " finished max_clique=" << overall_max_clique_size << " (" << total_branches << " branches in "
                << time_diff_str(begin, end) << ")\n";
    }

    INLINE bool enlarge_clique_size(const std::uint32_t thread_id, graph::size_t &max_clique_size, const graph::size_t depth) {
      bool found = false;

      {
        std::shared_lock clique_shared_lock(max_clique_mtx);
        if (depth > overall_max_clique_size) {
          clique_shared_lock.unlock();
          {
            std::scoped_lock clique_lock(max_clique_mtx);
            const graph::size_t old_overall_max_clique_size = overall_max_clique_size;
            overall_max_clique_size = std::max(depth, overall_max_clique_size);
            if (old_overall_max_clique_size != overall_max_clique_size) {
              owner_thread_id = thread_id;
            }

          }
          clique_shared_lock.lock();

          found = owner_thread_id == thread_id;

#ifdef LOG
          {
            std::scoped_lock print_lock(print_mtx);
            std::cerr << std::string((thread_id + 1), ' ') <<  thread_id << " enlarged clique to "
                      << overall_max_clique_size << "\n";
          }
#endif

        }
        max_clique_size = overall_max_clique_size;
      }

      return found;
    }

    void branch_heuristic(const std::uint32_t thread_id, const graph::key v,
                          std::vector<graph::key> &neighbours, std::vector<graph::key> &clique,
                          graph::size_t &max_clique_size, const graph::size_t upper_bound,
                          std::size_t &num_nodes, const graph::size_t depth = 1) {
      if (upper_bound_reached) {
        return;
      }

      num_nodes++;

      {
        std::shared_lock clique_shared_lock(max_clique_mtx);
        max_clique_size = overall_max_clique_size;
      }

      const graph::size_t prev_max_clique_size = max_clique_size;
      const graph::size_t next_depth = depth + 1;

      graph::key u = neighbours.back();
      neighbours.pop_back();

      std::vector<graph::key> new_neighbours = vertex_neighbourhood(u, neighbours);
      if (new_neighbours.empty() || next_depth == upper_bound) {
        if (enlarge_clique_size(thread_id, max_clique_size, next_depth)) {
          if (next_depth == upper_bound) {
            upper_bound_reached = true;
          }
          clique.clear();
          clique.reserve(next_depth);
          clique.emplace_back(u);
        }
      } else {
        std::vector<graph::colour> new_colours = sort_vertices_by_colour(new_neighbours);
        if (next_depth + new_colours.back() > max_clique_size) {
          branch_heuristic(thread_id, u, new_neighbours, clique, max_clique_size, upper_bound, num_nodes, next_depth);
        }
      }

      if (prev_max_clique_size < max_clique_size && owner_thread_id == thread_id) {
        clique.emplace_back(v);
      }
    }

    void branch_exact(const std::uint32_t thread_id, const graph::key v,
                      std::vector<graph::key> &neighbours, std::vector<graph::colour> &colours,
                      std::vector<graph::key> &clique, graph::size_t &max_clique_size,
                      const graph::size_t upper_bound,
                      std::size_t &num_nodes, const graph::size_t depth = 1) {
      num_nodes++;

      {
        std::shared_lock clique_shared_lock(max_clique_mtx);
        max_clique_size = overall_max_clique_size;
      }

      const graph::size_t next_depth = depth + 1;
      while (not neighbours.empty() && not upper_bound_reached && depth + colours.back() > max_clique_size) {
        const graph::size_t prev_max_clique_size = max_clique_size;

        graph::key u = neighbours.back();
        neighbours.pop_back();
        colours.pop_back();

        std::vector<graph::key> new_neighbours = vertex_neighbourhood(u, neighbours);
        if (new_neighbours.empty() || next_depth == upper_bound) {
          if (enlarge_clique_size(thread_id, max_clique_size, next_depth)) {
            if (next_depth == upper_bound) {
              upper_bound_reached = true;
            }
            clique.clear();
            clique.reserve(next_depth);
            clique.emplace_back(u);
          }
        } else {
          std::vector<graph::colour> new_colours = sort_vertices_by_colour(new_neighbours);
          if (next_depth + new_colours.back() > max_clique_size) {
            branch_exact(thread_id, u, new_neighbours, new_colours, clique, max_clique_size, upper_bound, num_nodes, next_depth);
          }
        }

        if (prev_max_clique_size < max_clique_size && owner_thread_id == thread_id) {
          clique.emplace_back(v);
        }

      }
    }

    INLINE std::vector<graph::key>
    vertex_neighbourhood(const graph::key v, const std::vector<graph::key> &neighbours) const {
#ifdef PROFILER
      Profiler::InstrumentationTimer timer_neighbourhood("Inducing The Neighbourhood");
#endif

      std::vector<graph::key> new_neighbours;
      new_neighbours.reserve(neighbours.size());

      for (graph::key u : neighbours) {
        if (adj_mtx[v][u]) {
          new_neighbours.emplace_back(u);
        }
      }

      return new_neighbours;
    }

    bool is_clique(const std::vector<graph::key> &clique) const {
      for (const auto &v : clique) {
        for (const auto &u : clique) {
          if (v != u) {
            if (auto N_u = neighbours_at(u); std::find(N_u.begin(), N_u.end(), v) == N_u.end()) {
              return false;
            }
          }
        }
      }
      return true;
    }
  };
}

int main(int ac, const char *av[]) {
  using namespace graph;

  std::cout << std::fixed << std::setprecision(3) << std::left;
  std::cerr << std::fixed << std::setprecision(3) << std::left;

#ifdef LOG
  const auto hold_time = 1s;
#endif

  std::cout << "Num Threads " << num_threads << "\n";

  try {
    fast_maxclique fmc;

#ifndef TEST_BOUNDS
    const std::size_t l = 0, r = -1u;
#else
    const std::size_t lb = 0, ub = fmc.expected_max_clique_size;
    for (std::size_t l = lb; l <= ub; ++l) {
      for (std::size_t r = l; r <= ub; ++r) {
#endif
        std::size_t num_turns = ac != 1 ? std::size_t(std::atol(av[1])) : -1u;
        while (num_turns--) {
//           std::cout << "Using Heuristic \n";
//           fmc.find_maxclique(fast_maxclique::algorithms::heuristic, l, r);
// #ifdef LOG
//           std::this_thread::sleep_for(hold_time);
// #endif

          std::cout << "Using Hybrid [" << l << ", " << r << "]\n";
          fmc.find_maxclique(fast_maxclique::algorithms::hybrid, l, r);
#ifdef LOG
          std::this_thread::sleep_for(hold_time);
#endif

          std::cout << "Using Exact [" << l << ", " << r << "]\n";
          fmc.find_maxclique(fast_maxclique::algorithms::exact, l, r);
#ifdef LOG
          std::this_thread::sleep_for(hold_time);
#endif

        }
#ifdef TEST_BOUNDS
      }
    }
#endif
  } catch (const std::exception &e) {
    std::cout << e.what() << "\n";
  }

  return 0;
}
