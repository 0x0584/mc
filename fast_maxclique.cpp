#include <assert.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
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

//#define NDEBUG

using namespace std::chrono_literals;

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

  INLINE end_of_scope_executor(Callable &&callback) : callback(std::move(callback)) {}
  INLINE ~end_of_scope_executor() { callback(); }

  end_of_scope_executor &operator=(const end_of_scope_executor &) = delete;
  end_of_scope_executor &operator=(end_of_scope_executor &&) = delete;

private:
  Callable callback;
};

namespace graph {
  using vertex = unsigned;
  using key = unsigned;
  using size_t = unsigned;
  using colour = unsigned;

  const std::uint32_t threads_per_core = 2;
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

    graph::size_t expected_max_clique_size, num_v, num_e;

    std::vector<graph::key> max_clique;
    std::vector<std::pair<graph::vertex, std::vector<graph::key>>> adj_lst; // Adjacency List
    std::vector<std::vector<bool>> adj_mtx; // Adjacency Matrix

#define neighbours_at(index)                    adj_lst[index].second

    fast_maxclique() {
      auto begin = std::chrono::high_resolution_clock::now();

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
      assert(tmp_adj_lst.size() == num_v);

      //  ------ sort using neighbour degrees
      //
      std::unordered_map<graph::vertex, graph::size_t> tmp_adj_lst_neighbours_degree;
      tmp_adj_lst_neighbours_degree.max_load_factor(0.5);
      tmp_adj_lst_neighbours_degree.reserve(num_v);
      for (const auto &[v, neighbours] : tmp_adj_lst) {
        for (const auto &u : neighbours) {
          tmp_adj_lst_neighbours_degree[v] += tmp_adj_lst[u].size();
        }
      }

      //  ------ sort using core numbers
      //
      // std::unordered_map<graph::vertex, size_t> tmp_adj_lst_cores;
      // tmp_adj_lst_cores.max_load_factor(0.5);
      // tmp_adj_lst_cores.reserve(num_v);
      // for (const auto &[v, neighbours] : tmp_adj_lst) {
      //   tmp_adj_lst_cores[v] = neighbours.size();
      // }

      // {
      //   size_t todo = tmp_adj_lst_cores.size();
      //   for (size_t level = 1; todo; ++level) {
      //     std::vector<graph::vertex> curr;
      //     curr.reserve(tmp_adj_lst_cores.size());

      //     for (const auto &[v, degree] : tmp_adj_lst_cores) {
      //       if (degree == level) {
      //         curr.emplace_back(v);
      //       }
      //     }

      //     for (std::vector<graph::vertex> next; not curr.empty(); curr = std::move(next)) {
      //       todo -= curr.size();
      //       next.reserve(curr.size());
      //       for (const auto &v : curr) {
      //         for (const auto &u : tmp_adj_lst[v]) {
      //           if (auto &k = tmp_adj_lst_cores[u]; k > level && --k == level) {
      //             next.emplace_back(u);
      //           }
      //         }
      //       }
      //     }
      //   }
      // }

      std::vector<std::pair<graph::vertex, std::set<graph::vertex>>> adj_lst_orig;
      adj_lst_orig.reserve(num_v);
      for (auto &v : tmp_adj_lst) {
        adj_lst_orig.emplace_back(std::make_pair(v.first, std::move(v.second)));
      }

      //  ------ sort using neighbour degrees
      //
      std::sort(adj_lst_orig.begin(), adj_lst_orig.end(),
                [&](const std::pair<graph::vertex, std::set<graph::vertex>> &u,
                    const std::pair<graph::vertex, std::set<graph::vertex>> &v) {
                  return (u.second.size() < v.second.size() ||
                          (u.second.size() == v.second.size() &&
                           tmp_adj_lst_neighbours_degree.at(u.first) < tmp_adj_lst_neighbours_degree.at(v.first)));
                });



      //  ------ sort using cores
      //
      // std::sort(adj_lst_orig.begin(), adj_lst_orig.end(),
      //           [&](const std::pair<graph::vertex, std::set<graph::vertex>> &u,
      //               const std::pair<graph::vertex, std::set<graph::vertex>> &v) {
      //             const auto &u_core = tmp_adj_lst_cores.at(u.first),
      //               &v_core = tmp_adj_lst_cores.at(v.first);
      //             const auto &u_degree = u.second.size(),
      //               &v_degree = v.second.size();
      //             const auto &u_neighbours_degree = tmp_adj_lst_neighbours_degree.at(u.first),
      //               &v_neighbours_degree = tmp_adj_lst_neighbours_degree.at(v.first);

      //             // return (u_degree < v_degree ||
      //             //         (u_degree == v_degree && u_neighbours_degree < v_neighbours_degree) ||
      //             //         (u_degree == v_degree && u_neighbours_degree == v_neighbours_degree &&
      //             //          u_cores < v_cores) ||
      //             //         (u_degree == v_degree && u_neighbours_degree == v_neighbours_degree &&
      //             //          u_cores == v_cores && u < v));
      //             return u_core < v_core || ( u_core == v_core && u_degree < v_degree) ||
      //               ( u_core == v_core && u_degree == v_degree &&
      //                 u_neighbours_degree < v_neighbours_degree);
      //           });


      //  ------ sort using vertex degrees
      //
      // std::sort(adj_lst_orig.begin(), adj_lst_orig.end(),
      //           [&](const std::pair<graph::vertex, std::set<graph::vertex>> &u,
      //               const std::pair<graph::vertex, std::set<graph::vertex>> &v) {
      //             return u.second.size() < v.second.size();
      //           });

      std::unordered_map<graph::vertex, graph::key> reversed_mapping;
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

      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Graph (vertices=" << num_v << " edges=" << num_e
                << ") with expected maxlique=" << expected_max_clique_size
                << " was read in " << std::chrono::duration<double>(end-begin).count() << "s\n";

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

      std::this_thread::sleep_for(3s);

    }

    std::vector<graph::vertex>
    find_maxclique(algorithms algo_type = algorithms::exact, graph::size_t upper_bound = -1u) {
      solution(algo_type, upper_bound);

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

      if (algo_type != algorithms::heuristic
          && expected_max_clique_size != clique.size()) {
        throw std::runtime_error(" NOT A MAX CLIQUE!!\n");
      }

      max_clique.clear();
      overall_max_clique_size = 0;

      return clique;
    }

  private:
#ifdef LOG
    mutable std::mutex print_mtx;
    void print(const std::string &s) const {
      std::scoped_lock print_lock(print_mtx);
      std::cerr << "> " << s << std::endl;
    }

    void print(const std::string &s, std::uint32_t thread_id) const {
      std::scoped_lock print_lock(print_mtx);
      std::cerr << (std::string((thread_id + 1), ' ') + std::to_string(thread_id))
                << " " << s << std::endl;
    }
#endif

    struct thread_wrapper {
      std::thread th;
      bool available = true;

      ~thread_wrapper() {
        if (th.joinable()) {
          th.join();
        }
      }

      thread_wrapper &operator=(std::thread &&rhs) {
        if (not available) {
          throw std::runtime_error("Thread is still running");
        }
        available = false;
        if (th.joinable()) {
          th.join();
        }
        th = std::move(rhs);
        return *this;
      }
    };

    INLINE graph::size_t num_colours_greedy(const std::vector<graph::key> &S) const {
      std::unordered_set<graph::colour> C;
      std::unordered_map<graph::key, graph::colour> L;

      L.max_load_factor(0.5);
      L.reserve(S.size());
      C.max_load_factor(0.5);
      C.reserve(S.size());

      for (const auto &v : S) {
        std::set<graph::colour> N_C;
        for (const auto &w : S) {
          for (const auto &u : neighbours_at(v)) {
            if (u == w) {
              N_C.emplace(L[u]);
              break;
            }
          }
        }
        N_C.erase(0);

        graph::colour colour = 1;
        for (const auto &c : N_C) {
          if (c != colour) {
            break;
          }
          colour++;
        }

        L[v] = colour;
        C.emplace(colour);
      }

      return C.size();
    }

    INLINE std::vector<graph::key> core_numbers(std::vector<graph::key> K) {
      graph::size_t todo = K.size();
      for (graph::size_t v = 0; v < K.size(); ++v) {
        if (K[v] == 0) {
          todo--;
        }
      }

      for (graph::size_t level = 1; todo; ++level) {
        std::vector<graph::key> curr;
        curr.reserve(K.size());

        for (graph::key v = 0; v < K.size(); ++v) {
          if (K[v] == level) {
            curr.emplace_back(v);
          }
        }

        for (std::vector<graph::key> next; not curr.empty(); curr = std::move(next)) {
          todo -= curr.size();
          next.reserve(curr.size());
          for (const auto &v : curr) {
            for (const auto &u : neighbours_at(v)) {
              if (K[u] > level && --K[u] == level) {
                next.emplace_back(u);
              }
            }
          }
        }
      }

      return K;
    }

    std::vector<graph::size_t> degrees() const {
      std::vector<graph::size_t> D;
      D.reserve(num_v);
      for (const auto &[_, neighbours] : adj_lst) {
        D.emplace_back(neighbours.size());
      }
      return D;
    }

    INLINE std::unordered_map<graph::key, graph::colour>
    sort_vertices_by_colour(std::vector<graph::key> &neighbours, const std::uint32_t thread_id) const {
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
        if (colour >= num_colours)
          num_colours = colour + 1;
      }
      // print("colouring " + std::to_string(neighbours.size()) + " vertices with "
      //       + std::to_string(num_colours) + " colours", thread_id);

      neighbours.clear();

      std::unordered_map<graph::key, graph::colour> colours;
      colours.max_load_factor(0.5);
      colours.reserve(neighbours.size());
      for (graph::colour colour = 0; colour < num_colours; ++colour) {
        for (const auto &v : colour_classes[colour]) {
          colours.emplace(v, colour + 1);
          neighbours.emplace_back(v);
        }
      }

      return colours;
    }

    mutable std::shared_mutex max_clique_mtx;
    mutable graph::size_t overall_max_clique_size = 0;
    void solution(algorithms algo_type, graph::size_t upper_bound) {
      if (algo_type == algorithms::hybrid) {
        solution(algorithms::heuristic, upper_bound);
        if (overall_max_clique_size >= upper_bound) {
          return;
        }
      }

      std::vector<graph::size_t> D = degrees();
      std::vector<graph::size_t> K = core_numbers(D);
      if (algo_type == algorithms::hybrid) {
        for (graph::key v = 0; v < num_v; ++v) {
          if (K[v] < overall_max_clique_size) {
            for (const auto &u : neighbours_at(v)) {
              if (D[u]) {
                D[u]--;
              }
            }
            K[v] = D[v] = 0;
          }
        }
        K = core_numbers(D);
      }

      ///////////////////////

      std::shared_mutex thread_mtx;
      std::condition_variable_any thread_is_available;
      std::vector<thread_wrapper> threads(num_threads);
      std::uint32_t owner_thread_id = -1u;

      std::mutex read_mtx;
      std::condition_variable reading_is_complete;
      bool reading = false;

      std::atomic_bool end_of_search = false;

      /////////////////////////

      auto begin = std::chrono::high_resolution_clock::now();
#ifndef LOG
      auto percent_begin = std::chrono::high_resolution_clock::now();
      graph::size_t old_max_clique_size = overall_max_clique_size;
#endif

      std::vector<graph::key> vertices;
      vertices.reserve(num_v);
      for (graph::key v = 0; v < num_v; ++v) {
        if (D[v] != 0) {
          vertices.emplace_back(v);
        }
      }

      std::unordered_map<graph::key, graph::colour> colours = sort_vertices_by_colour(vertices, 0);
      for (auto itr = vertices.rbegin(); itr != vertices.rend();) {
#ifdef LOG
        print("waiting... " + std::to_string(((std::distance(vertices.rbegin(), itr) + 1) * 100 / vertices.size())) + "%");
#endif

        //{
        std::unique_lock thread_lock(thread_mtx);
        thread_is_available.wait(thread_lock, [&] {
          return std::any_of(threads.begin(), threads.end(), [](const auto &th) {
            return th.available;
          });
        });
        //}

        graph::size_t current_max_clique_size;
        {
          std::shared_lock clique_lock(max_clique_mtx);
          current_max_clique_size = overall_max_clique_size;
        }
        if (current_max_clique_size >= upper_bound) {
          break;
        }

#ifdef LOG
        print("new thread(s) available! max_clique=" + std::to_string(current_max_clique_size));
#else
        if (old_max_clique_size != current_max_clique_size) {
          std::cerr << " clique_size=" << current_max_clique_size << "\n";
          old_max_clique_size = current_max_clique_size;
        }
#endif
        for (std::uint32_t thread_id = 0; thread_id < threads.size() && itr != vertices.rend(); ++thread_id) {
          if (threads[thread_id].available) {
            const graph::key v = *itr++;
            if (colours[v] <= current_max_clique_size) {
              end_of_search = true;
#ifdef LOG
              print("end of search max_clique=" + std::to_string(current_max_clique_size));
#endif
              break;
            }

            reading = true;
            threads[thread_id] = std::thread([&, this](const std::uint32_t thread_id, const graph::key v,
                                                       graph::size_t max_clique_size) {
              end_of_scope_executor set_thread_available([&threads, &thread_id,
                                                          &thread_mtx, &thread_is_available] {
                {
                  std::shared_lock thread_lock(thread_mtx);
                  threads[thread_id].available = true;
                }
                thread_is_available.notify_one();
              });

              std::unordered_map<graph::key, graph::size_t> K_;
              std::vector<graph::key> neighbours;
              K_.max_load_factor(0.5);
              K_.reserve(D[v]);
              neighbours.reserve(D[v]);

              graph::size_t max_core = K[v];
              for (const auto &u : neighbours_at(v)) {
                if (D[u] != 0) {
                  if (K[u] >= max_clique_size) {
                    K_.emplace(u, K[u]);
                    neighbours.emplace_back(u);
                    max_core = std::max(max_core, K[u]);
                  }
                  --D[u];
                }
              }
              D[v] = 0;

              {
                std::scoped_lock reading_lock(read_mtx);
                reading = false;
              }
              reading_is_complete.notify_one();

              if (max_core >= max_clique_size && neighbours.size() >= max_clique_size) {
#ifdef LOG
                print(std::string("branching on ") + std::to_string(adj_lst[v].first)
                      + " with " + std::to_string(neighbours.size()) + " neighbours"
                      //+ " Core=" +  std::to_string(max_core)
                      , thread_id);
                auto begin = std::chrono::high_resolution_clock::now();
#endif

                std::vector<graph::key> clique;
                if (algo_type == algorithms::heuristic) {
                  std::sort(neighbours.begin(), neighbours.end(),
                            [&](const std::size_t &u, const std::size_t &v) {
                              return K_.at(u) < K_.at(v);
                            });
                  branch_heuristic(v, neighbours, owner_thread_id, thread_id, max_clique_size, clique, K_, end_of_search);
                } else {
                  branch_exact(v, neighbours, owner_thread_id, thread_id, max_clique_size, clique, K_, end_of_search);
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

#ifdef LOG
                auto end = std::chrono::high_resolution_clock::now();
                {
                  std::scoped_lock print_lock(print_mtx);
                  std::cerr << std::string((thread_id + 1), ' ') <<  thread_id << " done with " << adj_lst[v].first
                            << " took " << std::chrono::duration<double>(end - begin).count() << "s\n";
                }
#endif

              }
            }, thread_id, v, current_max_clique_size);

            std::unique_lock reading_lock(read_mtx);
            reading_is_complete.wait(reading_lock, [&reading]() {
              return not reading;
            });

            K = core_numbers(D);
          }
        }

        if (end_of_search) {
          break;
        }

#ifndef LOG
        auto percent_end = std::chrono::high_resolution_clock::now();
        if (std::chrono::duration<double>(percent_end - percent_begin).count() >= 1.) {
          std::cerr << ((std::distance(vertices.rbegin(), itr) + 1) * 100 / vertices.size()) << "% "
                    << std::chrono::duration<double>(percent_end - begin).count() << "s\n";
        }
        percent_begin = std::chrono::high_resolution_clock::now();
#endif
      }

#ifdef LOG
      print("stopped making threads\n");
#endif
      threads.clear();

      auto end = std::chrono::high_resolution_clock::now();
#ifdef LOG
      std::cerr << "\n> all threads joined\n\n";
#endif
      std::cout << "Clique is found using " << algo_type << " in "
                << std::chrono::duration<double>(end - begin).count() << "s\n";
    }

    INLINE bool max_clique_found(std::vector<graph::key> &clique,
                                 std::uint32_t &owner_thread_id, const std::uint32_t thread_id,
                                 graph::size_t &max_clique_size, const graph::size_t depth) const {
      bool found = false;

      {
        std::shared_lock clique_shared_lock(max_clique_mtx);
        if (depth > overall_max_clique_size) {

#ifdef LOG
          auto begin = std::chrono::high_resolution_clock::now();
#endif

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
          auto end = std::chrono::high_resolution_clock::now();
          {
            std::scoped_lock print_lock(print_mtx);
            std::cerr << std::string((thread_id + 1), ' ') <<  thread_id << " clique update to "
                      << overall_max_clique_size << " took " << std::chrono::duration<double>(end - begin).count() << "s\n";
          }
#endif

        }
        max_clique_size = overall_max_clique_size;
      }

      if (found) {
        clique.clear();
        clique.reserve(depth);
      }

      return found;
    }

    void branch_heuristic(const graph::key v, const std::vector<graph::key> &neighbours,
                          std::uint32_t &owner_thread_id, const std::uint32_t thread_id,
                          graph::size_t &max_clique_size, std::vector<graph::key> &clique,
                          const std::unordered_map<graph::key, graph::size_t> &K_,
                          const std::atomic_bool &end_of_search,
                          const graph::size_t depth = 1) const {
      if (end_of_search && owner_thread_id != thread_id) {
        return;
      }
      if (neighbours.empty()) {
        if (max_clique_found(clique, owner_thread_id, thread_id, max_clique_size, depth)) {
          clique.emplace_back(v);
        }
        return;
      }

      const graph::key v_max = neighbours.back();
      const graph::size_t prev_max_clique_size = max_clique_size;
      branch_heuristic(v_max, vertex_neighbourhood(v_max, neighbours, max_clique_size, K_),
                       owner_thread_id, thread_id, max_clique_size, clique, K_, end_of_search, depth + 1);
      if (end_of_search && owner_thread_id != thread_id) {
        return;
      }
      if (prev_max_clique_size < max_clique_size && owner_thread_id == thread_id) {
        clique.emplace_back(v);
      }
    }

    void branch_exact(const graph::key v, std::vector<graph::key> &neighbours,
                      std::uint32_t &owner_thread_id, const std::uint32_t thread_id,
                      graph::size_t &max_clique_size, std::vector<graph::key> &clique,
                      const std::unordered_map<graph::key, graph::size_t> &K_,
                      const std::atomic_bool &end_of_search,
                      const graph::size_t depth = 1) const {
      if (end_of_search && owner_thread_id != thread_id) {
        return;
      }

      if (neighbours.empty()) {
        if (max_clique_found(clique, owner_thread_id, thread_id, max_clique_size, depth)) {
          clique.emplace_back(v);
        }
        return;
      }

      //if (thread_id == 1) {
      //print(" beginning at depth " + std::to_string(depth) + " with " + std::to_string(neighbours.size()), thread_id);
         //}
      std::unordered_map<graph::key, graph::colour> colours = sort_vertices_by_colour(neighbours, thread_id);
      while (not neighbours.empty() && depth + colours.at(neighbours.back()) > max_clique_size) {
        const graph::key u = neighbours.back();
        neighbours.pop_back();

        std::vector<graph::key> new_neighbours = vertex_neighbourhood(u, neighbours, max_clique_size, K_);
        const graph::size_t prev_max_clique_size = max_clique_size;
        branch_exact(u, new_neighbours, owner_thread_id, thread_id, max_clique_size, clique, K_, end_of_search, depth + 1);
        if (end_of_search && owner_thread_id != thread_id) {
          return;
        }
        if (prev_max_clique_size < max_clique_size && thread_id == owner_thread_id) {
          clique.emplace_back(v);
        }
      }
      // if (thread_id == 1) {
      // print(" done with depth " + std::to_string(depth), thread_id);
      // }
    }

    INLINE std::vector<graph::key>
    vertex_neighbourhood(const graph::key w, const std::vector<graph::key> &neighbours,
                         const graph::size_t max_clique_size,
                         const std::unordered_map<graph::key, graph::size_t> &K_) const {
      std::vector<graph::key> new_neighbours;
      new_neighbours.reserve(neighbours.size());

      for (const auto &u : neighbours) {
        if (std::binary_search(neighbours_at(w).begin(), neighbours_at(w).end(), u)
            && K_.at(u) >= max_clique_size) {
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

  void print_clique(const std::vector<graph::vertex> &) {
    // std::set<graph::vertex> m(mc.begin(), mc.end());
    // std::cout << "SIZE = " << m.size() << " { ";
    // for (const auto &e : m) {
    //   std::cout << e << " ";
    // }
    // std::cout << "}\n\n";
  }
}

int main(int ac, const char *av[]) {
  using namespace graph;

  std::cout << std::fixed << std::setprecision(3) << std::left;
  std::cerr << std::fixed << std::setprecision(3) << std::left;

  try {
    fast_maxclique fmc;
    //std::this_thread::sleep_for(3s);

    std::size_t num_turns = ac != 1 ? std::size_t(std::atol(av[1])) : -1u;
    while (num_turns--) {
      std::cout << "Using Heuristic \n";
      print_clique(fmc.find_maxclique(fast_maxclique::algorithms::heuristic));
      //std::this_thread::sleep_for(2s);

      std::cout << "Using Hybrid \n";
      print_clique(fmc.find_maxclique(fast_maxclique::algorithms::hybrid));
      //std::this_thread::sleep_for(2s);

      std::cout << "Using Exact \n";
      print_clique(fmc.find_maxclique(fast_maxclique::algorithms::exact));
      //std::this_thread::sleep_for(2s);
    }
  } catch (const std::exception &e) {
    std::cout << e.what() << "\n";
  }

  return 0;
}
