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

namespace graph {
  using vertex_t = unsigned long long;
  using edge_t = std::pair<vertex_t, vertex_t>;

  const std::uint32_t threads_per_core = 1;
  const std::uint32_t num_threads = std::thread::hardware_concurrency() * threads_per_core;

  std::size_t expected_max_clique_size = 0;

  INLINE std::vector<std::pair<vertex_t, std::set<vertex_t>>> read_adjacency_list() {
    std::size_t num_v, num_e;
    std::cin >> num_v >> num_e >> expected_max_clique_size;

    std::unordered_map<vertex_t, std::set<vertex_t>> tmp_adj_lst;
    tmp_adj_lst.max_load_factor(0.5);
    tmp_adj_lst.reserve(num_v);
    auto begin = std::chrono::high_resolution_clock::now();
    for (std::size_t n_edges = num_e; n_edges--;) {
      std::size_t u, v;
      std::cin >> u >> v;
      tmp_adj_lst[u].emplace(v);
      tmp_adj_lst[v].emplace(u);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Graph (vertices=" << num_v << " edges=" << num_e
              << ") with expected maxlique=" << graph::expected_max_clique_size
              << " was read in " << std::chrono::duration<double>(end-begin).count() << "s\n";
    assert(tmp_adj_lst.size() == num_v);

    std::vector<std::pair<vertex_t, std::set<vertex_t>>> adj_lst;
    adj_lst.reserve(tmp_adj_lst.size());
    for (auto &v : tmp_adj_lst) {
      adj_lst.emplace_back(std::make_pair(v.first, std::move(v.second)));
    }
    return adj_lst;
  }

  struct fast_maxclique {
#define neighbours_at(index)                    A[index].second

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

    fast_maxclique() {
      auto begin = std::chrono::high_resolution_clock::now();
      std::vector<std::pair<vertex_t, std::set<vertex_t>>> adj_lst = graph::read_adjacency_list();

      assert(not adj_lst.empty());
      std::sort(adj_lst.begin(), adj_lst.end(),
                [](const std::pair<vertex_t, std::set<vertex_t>> &u,
                   const std::pair<vertex_t, std::set<vertex_t>> &v) {
                  return u.second.size() < v.second.size();
                });

      std::unordered_map<vertex_t, std::size_t> reversed_mapping;
      reversed_mapping.reserve(adj_lst.size());
      A.reserve(adj_lst.size());
      for (const auto &v : adj_lst) {
        reversed_mapping.emplace(v.first, A.size());
        A.emplace_back(std::make_pair(v.first, std::vector<std::size_t>{}));
        A.back().second.reserve(v.second.size());
      }

      for (const auto &v : adj_lst) {
        const std::size_t &index = reversed_mapping[v.first];
        std::vector<std::size_t> &neighbours = A[index].second;

        for (const auto &u : v.second) {
          neighbours.emplace_back(reversed_mapping[u]);
        }

        std::sort(neighbours.begin(), neighbours.end(),
                  [&](const std::size_t &u, const std::size_t &v) {
                    return neighbours_at(u).capacity() < neighbours_at(v).capacity();
                  });
      }

      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Graph was prepared in "
                << std::chrono::duration<double>(end-begin).count() << "s\n\n";

#if 0//def LOG
      for (const auto &v : A) {
        std::cout << "Vertex: " << v.first << " { ";
        for (const auto &u : v.second) {
          std::cout << A[u].first << " ";
        }
        std::cout << "}\n";
      }
      std::cout << "\n";

      std::map<std::size_t, std::vector<vertex_t>> degrees;
      for (const auto &v : A) {
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

      //#ifndef LOG
      std::map<std::size_t, std::vector<vertex_t>> degrees;
      for (const auto &v : A) {
        degrees[v.second.size()].emplace_back(v.first);
      }
      for (const auto &deg : degrees) {
        std::cerr << "Degree " << deg.first << ": " << deg.second.size() << "\n";
      }
      std::cerr << "\n";
      //#endif

      std::this_thread::sleep_for(3s);

    }

    std::vector<vertex_t> find_maxclique(algorithms algo_type = algorithms::exact,
                                         std::size_t upper_bound = -1u) {
      solution(algo_type, upper_bound);

      std::vector<vertex_t> clique;
      clique.reserve(max_clique.size());
      for (const auto &v : max_clique) {
        clique.emplace_back(A[v].first);
      }

      std::set<vertex_t> m(clique.begin(), clique.end());
      std::cout << "SIZE = " << m.size() << " { ";
      for (const auto &e : m) {
        std::cout << e << " ";
      }
      std::cout << "}\n\n";

      if (not is_clique(max_clique)) {
        throw std::runtime_error(" NOT A CLIQUE!!\n");
      }

      if (algo_type != algorithms::heuristic
          && graph::expected_max_clique_size != clique.size()) {
        throw std::runtime_error(" NOT A MAX CLIQUE!!\n");
      }

      max_clique.clear();

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

    INLINE std::size_t num_colours_greedy(const std::vector<std::size_t> &S) const {
      std::unordered_set<std::size_t> C;
      std::unordered_map<std::size_t, std::size_t> L;

      L.max_load_factor(0.5);
      L.reserve(S.size());
      C.max_load_factor(0.5);
      C.reserve(S.size());

      for (const auto &v : S) {
        std::set<std::size_t> N_C;
        for (const auto &w : S) {
          for (const auto &u : neighbours_at(v)) {
            if (u == w) {
              N_C.emplace(L[u]);
              break;
            }
          }
        }
        N_C.erase(0);

        std::size_t colour = 1;
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

    INLINE std::vector<std::size_t> core_numbers(std::vector<std::size_t> K) {
      std::size_t todo = K.size();
      for (std::size_t v = 0; v < K.size(); ++v) {
        if (K[v] == 0) {
          todo--;
        }
      }

      for (std::size_t level = 1; todo; ++level) {
        std::vector<std::size_t> curr;
        curr.reserve(K.size());

        for (std::size_t v = 0; v < K.size(); ++v) {
          if (K[v] == level) {
            curr.emplace_back(v);
          }
        }

        for (std::vector<std::size_t> next; not curr.empty(); curr = std::move(next)) {
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

    std::vector<std::size_t> degrees() const {
      std::vector<std::size_t> D;
      D.reserve(A.size());
      for (const auto &[_, neighbours] : A) {
        D.emplace_back(neighbours.size());
      }
      return D;
    }

    mutable std::shared_mutex max_clique_mtx;
    mutable std::size_t overall_max_clique_size = 0;
    void solution(algorithms algo_type, std::size_t upper_bound) {
      if (algo_type == algorithms::hybrid) {
        solution(algorithms::heuristic, upper_bound);
        if (max_clique.size() >= upper_bound) {
          return;
        }
      }

      overall_max_clique_size = std::max(1ul, max_clique.size());

      std::vector<std::size_t> D = degrees();
      std::vector<std::size_t> K = core_numbers(D);
      if (algo_type == algorithms::hybrid) {
        for (std::size_t v = 0; v < A.size(); ++v) {
          if (K[v] < max_clique.size()) {
            for (const auto &u : neighbours_at(v)) {
              if (D[u]) {
                D[u]--;
              }
            }
            D[v] = 0;
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

      /////////////////////////

      auto begin = std::chrono::high_resolution_clock::now();
#ifndef LOG
      auto percent_begin = std::chrono::high_resolution_clock::now();
      std::size_t old_max_clique_size = 0;
#endif

      for (std::size_t v = 0; v < A.size();) {
#ifdef LOG
        print("waiting... " + std::to_string(((v+1)*100 / A.size())) + "%");
#endif

        //{
          std::unique_lock thread_lock(thread_mtx);
          thread_is_available.wait(thread_lock, [&] {
            return std::any_of(threads.begin(), threads.end(), [](const auto &th) {
              return th.available;
            });
          });
          //}

        std::size_t current_max_clique_size;
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

        for (std::size_t thread_id = 0; thread_id < threads.size() && v < A.size(); ++thread_id) {
          if (threads[thread_id].available) {
            reading = true;
            threads[thread_id] = std::thread([&, this](std::uint32_t thread_id, std::size_t v,
                                                       std::size_t max_clique_size) {
              std::unordered_map<std::size_t, std::size_t> K_;
              std::vector<std::size_t> neighbours;
              K_.max_load_factor(0.5);
              K_.reserve(D[v]);
              neighbours.reserve(D[v]);
              std::size_t max_core = K[v];
              for (const auto &u : neighbours_at(v)) {
                if (K[u] >= max_clique_size) {
                  K_[u] = K[u];
                  neighbours.emplace_back(u);
                  max_core = std::max(max_core, K[u]);
                  assert(D[u] > 0);
                  --D[u];
                } else if (D[u] != 0) {
                  D[u] = 0;
                }
              }
              D[v] = 0;

              {
                std::scoped_lock reading_lock(read_mtx);
                reading = false;
              }
              reading_is_complete.notify_one();

              if (std::min(max_core, neighbours.size()) >= max_clique_size) {
#ifdef LOG
                print(std::string("branching on ") + std::to_string(A[v].first)
                      //+ " Core=" +  std::to_string(max_core)
                      , thread_id);
                auto begin = std::chrono::high_resolution_clock::now();
#endif

                std::sort(neighbours.begin(), neighbours.end(),
                          [&](const std::size_t &u, const std::size_t &v) {
                            return K_.at(u) < K_.at(v);
                          });

                std::vector<std::size_t> clique;

                if (algo_type == algorithms::heuristic) {
                  branch_heuristic(v, neighbours, owner_thread_id, thread_id, max_clique_size, clique, K_);
                } else {
                  branch_exact(v, neighbours, owner_thread_id, thread_id, max_clique_size, clique, K_);
                }

                if (std::scoped_lock clique_lock(max_clique_mtx);
                    clique.size() > max_clique.size() && thread_id == owner_thread_id) {
#ifdef LOG
                  std::ostringstream oss;
                  oss << "found clique for " << A[v].first << " of size=" << clique.size() << " { ";
                  for (const auto &v : clique) {
                    oss << A[v].first << " ";
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
                  std::cerr << std::string((thread_id + 1), ' ') <<  thread_id << " done with " << A[v].first
                            << " took " << std::chrono::duration<double>(end - begin).count() << "s\n";
                }
#endif

              }

              {
                std::shared_lock thread_lock(thread_mtx);
                threads[thread_id].available = true;
              }
              thread_is_available.notify_one();
            }, thread_id, v, current_max_clique_size);

            std::unique_lock reading_lock(read_mtx);
            reading_is_complete.wait(reading_lock, [&reading]() {
              return not reading;
            });

            K = core_numbers(D);
            v++;
          }
        }

#ifndef LOG
        auto percent_end = std::chrono::high_resolution_clock::now();
        if (v == num_threads || std::chrono::duration<double>(percent_end - percent_begin).count() >= 1.) {
          std::cerr << ((v+1)*100 / A.size()) << "% "
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

    INLINE bool max_clique_found(std::vector<std::size_t> &clique,
                                 std::uint32_t &owner_thread_id, const std::uint32_t thread_id,
                                 std::size_t &max_clique_size, const std::size_t depth) const {
      bool found;
      {
        std::shared_lock clique_shared_lock(max_clique_mtx);
        if ((found = depth > overall_max_clique_size)) {

#ifdef LOG
          auto begin = std::chrono::high_resolution_clock::now();
#endif

          clique_shared_lock.unlock();
          {
            std::scoped_lock clique_lock(max_clique_mtx);
            const std::size_t old_overall_max_clique_size = overall_max_clique_size;
            overall_max_clique_size = std::max(depth, overall_max_clique_size);
            if (old_overall_max_clique_size != overall_max_clique_size) {
              owner_thread_id = thread_id;
            }

          }
          clique_shared_lock.lock();
          max_clique_size = overall_max_clique_size;

#ifdef LOG
          auto end = std::chrono::high_resolution_clock::now();
          {
            std::scoped_lock print_lock(print_mtx);
            std::cerr << std::string((thread_id + 1), ' ') <<  thread_id << " clique update to "
                      << depth << " took " << std::chrono::duration<double>(end - begin).count() << "s\n";
          }
#endif

        }
      }

      if (found) {
        clique.clear();
        clique.reserve(depth);
      }
      return found;
    }

    void branch_heuristic(const std::size_t v, const std::vector<std::size_t> &neighbours,
                          std::uint32_t &owner_thread_id, const std::uint32_t thread_id,
                          std::size_t &max_clique_size, std::vector<std::size_t> &clique,
                          const std::unordered_map<std::size_t, std::size_t> &K_,
                          const std::size_t depth = 1) const {
      if (neighbours.empty()) {
        if (max_clique_found(clique, owner_thread_id, thread_id,
                             max_clique_size, depth)) {
          clique.emplace_back(v);
        }
        return;
      }

      const std::size_t v_max = neighbours.back();
      const std::size_t prev_max_clique_size = max_clique_size;
      branch_heuristic(v_max, vertex_neighbourhood(v_max, neighbours, max_clique_size, K_),
                       owner_thread_id, thread_id, max_clique_size, clique, K_, depth + 1);
      if (prev_max_clique_size < max_clique_size && owner_thread_id == thread_id) {
        clique.emplace_back(v);
      }
    }

    void branch_exact(const std::size_t v, std::vector<std::size_t> &neighbours,
                      std::uint32_t &owner_thread_id, const std::uint32_t thread_id,
                      std::size_t &max_clique_size, std::vector<std::size_t> &clique,
                      const std::unordered_map<std::size_t, std::size_t> &K_,
                      const std::size_t depth = 1) const {
      if (neighbours.empty()) {
        if (max_clique_found(clique, owner_thread_id, thread_id,
                             max_clique_size, depth)) {
          clique.emplace_back(v);
        }
        return;
      }

      while (not neighbours.empty() && depth + neighbours.size() > max_clique_size) {
        const std::size_t u = neighbours.back();
        neighbours.pop_back();

        std::vector<std::size_t> new_neighbours = vertex_neighbourhood(u, neighbours, max_clique_size, K_);
        if (depth + new_neighbours.size() >= max_clique_size) {
          const std::size_t prev_max_clique_size = max_clique_size;
          branch_exact(u, new_neighbours, owner_thread_id, thread_id, max_clique_size, clique, K_, depth + 1);
          if (prev_max_clique_size < max_clique_size && thread_id == owner_thread_id) {
            clique.emplace_back(v);
          }
        }
      }
    }

    INLINE std::vector<std::size_t>
    vertex_neighbourhood(std::size_t w, const std::vector<std::size_t> &neighbours,
                         const std::size_t max_clique_size,
                         const std::unordered_map<std::size_t, std::size_t> &K_) const {
      std::vector<std::size_t> new_neighbours;
      new_neighbours.reserve(neighbours.size());

      for (const auto &u : neighbours) {
        for (const auto &v : neighbours_at(w)) {
          if (v == u  && K_.at(v) >= max_clique_size) {
            new_neighbours.emplace_back(u);
          }
        }
      }

      return new_neighbours;
    }

    bool is_clique(const std::vector<std::size_t> &clique) const {
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

    std::vector<std::size_t> max_clique;
    std::vector<std::pair<vertex_t, std::vector<std::size_t>>> A; // Adjacency List
  };

  void print_clique(const std::vector<vertex_t> &) {
    // std::set<vertex_t> m(mc.begin(), mc.end());
    // std::cout << "SIZE = " << m.size() << " { ";
    // for (const auto &e : m) {
    //   std::cout << e << " ";
    // }
    // std::cout << "}\n\n";
  }
}

int main(int ac, const char *av[]) {
  using namespace graph;

  std::cout << std::fixed << std::setw(6) << std::setprecision(3) << std::left;
  std::cerr << std::fixed << std::setw(6) << std::setprecision(3) << std::left;

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
