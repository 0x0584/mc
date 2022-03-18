#include <assert.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <execution>
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

#include "profiler.h"

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

  const std::uint32_t threads_per_core = 1;
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

      std::unordered_map<graph::vertex, graph::size_t> neighbours_degree;
      neighbours_degree.max_load_factor(0.5);
      neighbours_degree.reserve(num_v);
      for (const auto &[v, neighbours] : tmp_adj_lst) {
        for (const auto &u : neighbours) {
          neighbours_degree[v] += tmp_adj_lst[u].size();
        }
      }

      std::vector<std::pair<graph::vertex, std::set<graph::vertex>>> orig_adj_lst;
      orig_adj_lst.reserve(num_v);
      for (auto &v : tmp_adj_lst) {
        orig_adj_lst.emplace_back(std::make_pair(v.first, std::move(v.second)));
      }

      std::sort(std::execution::par_unseq, orig_adj_lst.begin(), orig_adj_lst.end(),
                [&](const std::pair<graph::vertex, std::set<graph::vertex>> &v,
                    const std::pair<graph::vertex, std::set<graph::vertex>> &u) {
                  return (u.second.size() < v.second.size() ||
                          (u.second.size() == v.second.size() &&
                           neighbours_degree.at(u.first) < neighbours_degree.at(v.first)));
                });

      std::unordered_map<graph::vertex, graph::key> reversed_mapping;
      reversed_mapping.reserve(num_v);
      adj_lst.reserve(num_v);
      adj_mtx.reserve(num_v);
      for (const auto &v : orig_adj_lst) {
        reversed_mapping.emplace(v.first, adj_lst.size());
        adj_mtx.emplace_back(std::vector<bool>(num_v));
        adj_lst.emplace_back(std::make_pair(v.first, std::vector<graph::key>{}));
        adj_lst.back().second.reserve(v.second.size());
      }

      for (const auto &[v, neighbours_orig] : orig_adj_lst) {
        const graph::key &v_key = reversed_mapping.at(v);
        std::vector<bool> &v_adj_mtx = adj_mtx[v_key];
        std::vector<graph::key> &neighbours = neighbours_at(v_key);
        for (const auto &u : neighbours_orig) {
          const graph::key &u_key = reversed_mapping.at(u);
          v_adj_mtx[u_key] = true;
          neighbours.emplace_back(u_key);
        }
        std::sort(std::execution::par_unseq, neighbours.begin(), neighbours.end());
      }

      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Graph (vertices=" << num_v << " edges=" << num_e
                << ") with expected maxlique=" << expected_max_clique_size
                << " was read in " << std::chrono::duration<double>(end-begin).count() << "s\n";
#ifdef LOG
      std::this_thread::sleep_for(3s);
#endif
    }

    std::vector<graph::vertex>
    find_maxclique(algorithms algo_type = algorithms::exact, graph::size_t upper_bound = -1u) {
      Profiler::beginSession("Max Clique", "solution.json");
      solution(algo_type, upper_bound);
      Profiler::endSession();

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
      end_of_search = false;

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

      INLINE thread_wrapper &operator=(std::thread &&rhs) {
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

    INLINE std::vector<graph::size_t>
    core_numbers_neighbours(std::vector<graph::key> &neighbours,
                            const graph::size_t max_clique_size) const {

      if (neighbours.empty() || neighbours.size() < max_clique_size) {
        return {};
      }

      Profiler::InstrumentationTimer timer_core("Computing cores");

      std::vector<graph::size_t> K;
      K.resize(neighbours.size());

      std::unordered_map<graph::key, std::set<graph::key>> local_adj_lst;
      local_adj_lst.max_load_factor(0.5);
      local_adj_lst.reserve(neighbours.size());
      for (std::size_t i = 0; i < neighbours.size(); ++i) {
        for (std::size_t j = 0; j < neighbours.size(); ++j) {
          if (i != j) {
            if (adj_mtx[neighbours[i]][neighbours[j]]) {
              local_adj_lst[i].emplace(j);
              K[i]++;
            }
          }
        }
      }

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
            for (const auto &u : local_adj_lst.at(v)) {
              if (K[u] > level && --K[u] == level) {
                next.emplace_back(u);
              }
            }
          }
        }
      }

      std::unordered_map<graph::key, graph::size_t> neighbours_mapping;
      neighbours_mapping.max_load_factor(0.5);
      neighbours_mapping.reserve(neighbours.size());

      for (graph::size_t v = 0, size = neighbours.size(); v < size; ++v) {
        if (K.back() < max_clique_size) {
          while (not K.empty() && K.back() < max_clique_size) {
            neighbours.pop_back();
            K.pop_back();
            size--;
          }

          if (K.empty()) {
            return {};
          }

          assert(&neighbours[v] != &neighbours.back());

          std::swap(neighbours[v], neighbours.back());
          neighbours.pop_back();

          neighbours_mapping.emplace(neighbours[v], v);

          std::swap(K[v], K.back());
          K.pop_back();

          size--;

        }
      }


      std::sort(neighbours.begin(), neighbours.end(),
                [&K, &neighbours_mapping](const graph::key &u, const graph::key &v) {
                  return K[neighbours_mapping[u]] < K[neighbours_mapping[v]];
                });
      std::sort(K.begin(), K.end());

      // {
      //   std::scoped_lock print_lock(print_mtx);
      //   for (std::size_t i = 0; i < neighbours.size(); ++i) {
      //     std::cerr << "vertex " << neighbours[i] << " core=" << K[i] << "\n";
      //   }
      // }

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

    std::shared_mutex max_clique_mtx;
    graph::size_t overall_max_clique_size = 0;
    std::atomic_bool end_of_search = false;
    std::uint32_t owner_thread_id = -1u;

    void solution(algorithms algo_type, graph::size_t upper_bound) {
      if (algo_type == algorithms::hybrid) {
        solution(algorithms::heuristic, upper_bound);
        if (overall_max_clique_size >= upper_bound) {
          return;
        }
      }

      std::vector<graph::size_t> D = degrees();

      ///////////////////////

      std::shared_mutex thread_mtx;
      std::condition_variable_any thread_is_available;
      std::vector<thread_wrapper> threads(num_threads);

      std::mutex read_mtx;
      std::condition_variable reading_is_complete;
      bool reading = false;

      /////////////////////////

      auto begin = std::chrono::high_resolution_clock::now();
#ifndef LOG
      auto percent_begin = std::chrono::high_resolution_clock::now();
      graph::size_t old_max_clique_size = overall_max_clique_size;
#endif

      for (graph::key v = 0; v < num_v; ++v) {
#ifdef LOG
        print("waiting... " + std::to_string(((v + 1) * 100 / num_v)) + "%");
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
        for (std::uint32_t thread_id = 0; thread_id < threads.size() && v < adj_lst.size(); ++thread_id) {
          if (threads[thread_id].available) {
//             if (colours[v] <= current_max_clique_size) {
//               end_of_search = true;
// #ifdef LOG
//               print("end of search max_clique=" + std::to_string(current_max_clique_size));
// #endif
//               break;
//             }

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

              std::vector<graph::key> neighbours;
              neighbours.reserve(D[v]);
              for (const auto &u : neighbours_at(v)) {
                if (D[u] != 0) { // remove if < max_clique_size
                  --D[u];
                  neighbours.emplace_back(u);
                }
              }
              D[v] = 0;

              {
                std::scoped_lock reading_lock(read_mtx);
                reading = false;
              }
              reading_is_complete.notify_one();

              std::vector<graph::key> K = core_numbers_neighbours(neighbours, max_clique_size);
              if (K.empty()) {
                return;
              }

#ifdef LOG
              print(std::string("branching on ") + std::to_string(adj_lst[v].first)
                    + " with " + std::to_string(neighbours.size()) + " neighbours"
                    + " core=" +  std::to_string(K.back())
                    , thread_id);
              auto begin = std::chrono::high_resolution_clock::now();
#endif

              std::vector<graph::key> clique;
              if (algo_type == algorithms::heuristic) {
                branch_heuristic(thread_id, v, neighbours, clique, max_clique_size);
              } else {
                branch_exact(thread_id, v, neighbours, K, clique, max_clique_size);
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

            }, thread_id, v, current_max_clique_size);

            std::unique_lock reading_lock(read_mtx);
            reading_is_complete.wait(reading_lock, [&reading]() {
              return not reading;
            });
          }
        }

        if (end_of_search) {
          break;
        }

#ifndef LOG
        auto percent_end = std::chrono::high_resolution_clock::now();
        if (std::chrono::duration<double>(percent_end - percent_begin).count() >= 1.) {
          std::cerr << std::to_string(((v + 1) * 100 / num_v)) << "%"
                    << std::chrono::duration<double>(percent_end - begin).count() << "s\n";
        }
        percent_begin = std::chrono::high_resolution_clock::now();
#endif
      }

#ifdef LOG
      print("stopped making threads");
#endif
      threads.clear();

      auto end = std::chrono::high_resolution_clock::now();
#ifdef LOG
      std::cerr << "> all threads joined\n\n";
#endif
      std::cout << "Clique is found using " << algo_type << " in "
                << std::chrono::duration<double>(end - begin).count() << "s\n";
    }

    INLINE bool found_extended_clique(const std::uint32_t thread_id, graph::size_t &max_clique_size,
                                      const graph::size_t depth) {
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

      return found && not (end_of_search && owner_thread_id != thread_id);
    }

    void branch_heuristic(const std::uint32_t thread_id, const graph::key v,
                          std::vector<graph::key> &neighbours, std::vector<graph::key> &clique,
                          graph::size_t &max_clique_size, graph::size_t depth = 1) {
      if (end_of_search && owner_thread_id != thread_id) {
        return;
      }

      {
        std::shared_lock clique_shared_lock(max_clique_mtx);
        max_clique_size = std::max(overall_max_clique_size, max_clique_size);
      }

      const graph::size_t prev_max_clique_size = max_clique_size;

      graph::key u = neighbours.back();
      neighbours.pop_back();

      std::vector<graph::key> new_neighbours = neighbours;
      std::vector<graph::size_t> K = core_numbers_neighbours(new_neighbours, max_clique_size);

      if (K.empty()) {
        if (found_extended_clique(thread_id, max_clique_size, depth + 1)) {
          clique.clear();
          clique.reserve(depth + 1);
          clique.emplace_back(u);
        }
      } else if (depth + K.back() >= max_clique_size) {
        branch_heuristic(thread_id, u, new_neighbours, clique, max_clique_size, depth + 1);
      }

      if (prev_max_clique_size < max_clique_size && owner_thread_id == thread_id) {
        clique.emplace_back(v);
      }
    }

    void branch_exact(const std::uint32_t thread_id, const graph::key v,
                      std::vector<graph::key> &neighbours, std::vector<graph::key> &K,
                      std::vector<graph::key> &clique, graph::size_t &max_clique_size,
                      graph::size_t depth = 1) {
      {
        std::shared_lock clique_shared_lock(max_clique_mtx);
        max_clique_size = std::max(overall_max_clique_size, max_clique_size);
      }

#ifdef LOG
      auto begin = std::chrono::high_resolution_clock::now();
#endif

      while (not (end_of_search && owner_thread_id != thread_id) &&
             not neighbours.empty() && depth + K.back() >= max_clique_size) {
        const graph::size_t prev_max_clique_size = max_clique_size;

        graph::key u = neighbours.back();
        neighbours.pop_back();
        K.pop_back();

        std::vector<graph::key> new_neighbours;
        new_neighbours.reserve(neighbours.size());
        for (const auto &w : neighbours) {
          if (adj_mtx[w][u]) {
            new_neighbours.emplace_back(w);
          }
        }

        std::vector<graph::size_t> new_K = core_numbers_neighbours(new_neighbours, max_clique_size - depth);

        if (new_K.empty()) {
          if (found_extended_clique(thread_id, max_clique_size, depth + 1)) {
            clique.clear();
            clique.reserve(depth + 1);
            clique.emplace_back(u);
          }
        } else {
          branch_exact(thread_id, u, new_neighbours, new_K, clique, max_clique_size, depth + 1);
        }

        if (prev_max_clique_size < max_clique_size && thread_id == owner_thread_id) {
          clique.emplace_back(v);
        }

      }

#ifdef LOG
      auto end = std::chrono::high_resolution_clock::now();
      const double diff = std::chrono::duration<double>(end - begin).count();
      if (diff >= 1.) {
        std::scoped_lock print_lock(print_mtx);
        std::cerr << std::string((thread_id + 1), ' ') <<  thread_id << " end of branching on " << adj_lst[v].first
                  << " at depth " << depth << " took " << diff << "s\n";
      }
#endif

    }

    INLINE std::vector<graph::key>
    vertex_neighbourhood(const graph::key w, const std::vector<graph::key> &neighbours,
                         const graph::size_t max_clique_size,
                         //const std::unordered_map<graph::key, graph::size_t> &K_
                         const std::unordered_map<graph::key, graph::colour> &colours
                         ) {
      std::vector<graph::key> new_neighbours;
      new_neighbours.reserve(neighbours.size());

      for (const auto &u : neighbours) {
        // for (const auto &w neighbours_at(w)) {
          if (
              std::binary_search(neighbours_at(w).begin(), neighbours_at(w).end(), u)

              // && colours.at(u) >= max_clique_size
            ) {
            new_neighbours.emplace_back(u);
          }
        // }
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
      // std::cout << "Using Heuristic \n";
      // print_clique(fmc.find_maxclique(fast_maxclique::algorithms::heuristic));
      //std::this_thread::sleep_for(2s);

      // std::cout << "Using Hybrid \n";
      // print_clique(fmc.find_maxclique(fast_maxclique::algorithms::hybrid));
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
