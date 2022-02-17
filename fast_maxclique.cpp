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

using namespace std::chrono_literals;

#define LITERAL(expr) #expr

namespace graph {
  using vertex_t = unsigned long long;
  using edge_t = std::pair<vertex_t, vertex_t>;

  std::size_t expected_max_clique_size = 0;

  std::vector<std::pair<vertex_t, std::set<vertex_t>>> read_adjacency_list() {
    std::size_t num_v, num_e;
    std::cin >> num_v >> num_e >> expected_max_clique_size;

    std::map<vertex_t, std::set<vertex_t>> tmp_adj_lst;
    while (num_e--) {
      std::size_t u, v;
      std::cin >> u >> v;
      tmp_adj_lst[u].emplace(v);
      tmp_adj_lst[v].emplace(u);
    }

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
                  return u.second.size() > v.second.size();
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
                  [this](const std::size_t &u, const std::size_t &v) {
                    return neighbours_at(u).capacity() < neighbours_at(v).capacity();
                  });
      }

      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Graph is read in "
                << std::chrono::duration<double>(end-begin).count() << "s"
                << " with expected Max Clique " << graph::expected_max_clique_size << "\n\n";

#ifdef LOG
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
        // std::cout << "Vertex " << v.first << " has degree " << v.second.size() << "\n";
        // std::size_t prev_size = A[v.second.front()].second.size();
        // for (const auto &u : v.second) {
        //   if (prev_size < A[u].second.size()) {
        //     std::cout << v.first << " --- " << A[u].first << "\n";
        //     throw std::runtime_error(" SORTING IS NOT CORRECT ");
        //   }
        //   prev_size = A[u].second.size();
        // }
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

      //exit(1);
    }

    // branch iterative

    // clique.emplace_back(v);
    // while (not core_neighbours.empty()) {
    //   std::size_t v_max = max_degree_vertex(core_neighbours);
    //   clique.emplace_back(v_max);

    //   const std::size_t max_clique_size = std::max(clique.size(), max_clique.size());
    //   std::vector<std::size_t> new_neighbours;
    //   new_neighbours.reserve(neighbours_at(v_max).size());
    //   for (const auto &v : neighbours_at(v_max)) {
    //     if (not pruned[v] && neighbours_at(v).size() >= max_clique_size) {
    //       for (const auto &u : core_neighbours) {
    //         if (v == u) {
    //           new_neighbours.emplace_back(u);
    //         }
    //       }
    //     }
    //   }

    //   core_neighbours = std::move(new_neighbours);
    // }

    std::vector<vertex_t> find_maxclique(algorithms algo_type = algorithms::exact,
                                         std::size_t upper_bound = -1) {
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
    mutable std::mutex print_mtx;
    void print(const std::string &s, std::size_t thread_id = -1) const {
#ifdef LOG
      std::scoped_lock print_lock(print_mtx);
      std::cerr << (thread_id == size_t(-1)
                    ? ">"
                    : (std::string((thread_id + 1), ' ') + std::to_string(thread_id)))
                << " " << s << std::endl;
      //std::this_thread::sleep_for(0.05s);
#endif
    }

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

    void solution(algorithms algo_type, std::size_t upper_bound) {
      std::vector<bool> pruned(A.size());

      if (algo_type == algorithms::hybrid) {
        solution(algorithms::heuristic, upper_bound);
        if (max_clique.size() >= upper_bound) {
          return;
        }
      }

      ///////////////////////

      const std::uint32_t num_threads = std::thread::hardware_concurrency() * 5;

      std::shared_mutex thread_mtx;
      std::condition_variable_any thread_is_available;
      std::vector<thread_wrapper> threads(num_threads);
      std::vector<std::size_t> available_threads;
      available_threads.reserve(num_threads);

      std::mutex max_clique_mtx;
      std::size_t overall_max_clique_size = max_clique.size();

      std::mutex prune_mtx;
      std::condition_variable apply_prune;
      std::vector<std::size_t> to_prune;
      to_prune.reserve(num_threads);

      const auto valid_vertex = [&, this](std::size_t v) {
        return v < A.size() && overall_max_clique_size < neighbours_at(v).size();
      };

      /////////////////////////

      auto begin = std::chrono::high_resolution_clock::now();
      for (std::size_t v = 0; valid_vertex(v); available_threads.clear()) {
        print("waiting...");

        std::unique_lock thread_lock(thread_mtx);
        thread_is_available.wait(thread_lock, [&] {
          bool flag = false;

          for (std::uint32_t thread_id = 0; thread_id < threads.size()
                 && valid_vertex(v); ++thread_id) {
            if (threads[thread_id].available) {
              flag = true;
              available_threads.emplace_back(thread_id);
              to_prune.emplace_back(v++);
              print( "available", thread_id);
            }
          }

          return flag;
        });

        if (std::scoped_lock clique_lock(max_clique_mtx);
            (overall_max_clique_size = max_clique.size()) >= upper_bound) {
          break;
        }

        print(std::to_string(available_threads.size()) + " threads available!");

        std::uint32_t vertices_read = 0;
        auto current = to_prune.rbegin();
        for (const auto &thread_id : available_threads) {
          std::size_t v = *current++;

          print("branching on " + std::to_string(A[v].first), thread_id);
          threads[thread_id] =
            std::thread([&, this](std::uint32_t thread_id, std::size_t v,
                                  std::size_t max_clique_size,
                                  std::vector<std::size_t> to_prune) {
              std::vector<std::size_t> neighbours;
              neighbours.reserve(neighbours_at(v).size());
              for (const auto &u : neighbours_at(v)) {
                if (not pruned[u] && neighbours_at(u).size() > overall_max_clique_size
                    // && neighbours_at(u).size() >= neighbours_at(v).size()
                    && std::find(to_prune.begin(), to_prune.end(), u) == to_prune.end()) {
                  neighbours.emplace_back(u);
                }
              }

              {
                std::scoped_lock prune_lock(prune_mtx);
                if (++vertices_read == available_threads.size()) {
                  apply_prune.notify_one();
                }
              }

              std::vector<std::size_t> clique;
              const std::size_t prev_max_clique_size = max_clique_size;

              auto begin = std::chrono::high_resolution_clock::now();
              if (algo_type == algorithms::heuristic) {
                branch_heuristic(v, neighbours, max_clique_size, clique);
              } else {
                branch_exact(v, neighbours, max_clique_size, clique, thread_id);
              }

              if (max_clique_size > prev_max_clique_size) {
                std::scoped_lock clique_lock(max_clique_mtx);
#ifdef LOG
                std::ostringstream oss;
                oss << (max_clique_size > max_clique.size() ? "found" : "droped")
                    << " clique for " << A[v].first << " of size=" << max_clique_size << " { ";
                for (const auto &v : clique) {
                  oss << A[v].first << " ";
                }
                oss << "}";
                print(oss.str(), thread_id);
#endif

                if (max_clique_size > max_clique.size()) {
                  max_clique = std::move(clique);
                }

              }
#ifdef LOG
              auto end = std::chrono::high_resolution_clock::now();
              {
                std::scoped_lock print_lock(print_mtx);
                std::cerr << std::string((thread_id + 1), ' ') <<  thread_id << " done with " << A[v].first
                          << " took " << std::chrono::duration<double>(end - begin).count() << "s\n";
              }
#endif

              {
                std::shared_lock thread_lock(thread_mtx);
                threads[thread_id].available = true;
              }
              thread_is_available.notify_one();
            }, thread_id, v, overall_max_clique_size, std::vector<std::size_t>(current, to_prune.rend()));
        }

        {
          std::unique_lock prune_lock(prune_mtx);
          apply_prune.wait(prune_lock, [&]() {
            return vertices_read == available_threads.size();
          });
          for (const auto &v : to_prune) {
            pruned[v] = true;
          }
          to_prune.clear();
        }

      }

      print("stopped making threads\n");

      threads.clear();

      auto end = std::chrono::high_resolution_clock::now();
#ifdef LOG
      std::cerr << "\n> all threads joined\n\n";
#endif
      std::cout << "Clique is found using " << algo_type << " in "
                << std::chrono::duration<double>(end - begin).count() << "s\n";
    }

    void branch_heuristic(std::size_t v, const std::vector<std::size_t> &neighbours,
                          std::size_t &max_clique_size, std::vector<std::size_t> &clique,
                          std::size_t depth = 1) const {
      if (neighbours.empty()) {
        if (depth > max_clique_size) {
          max_clique_size = depth;
          clique.reserve(depth);
          clique.emplace_back(v);
        }
        return;
      }

      // std::size_t v_max = neighbours.front();
      // for (const auto &v : neighbours) {
      //   if (neighbours_at(v_max).size() < neighbours_at(v).size()) {
      //     v_max = v;
      //   }
      // }

      // assert(neighbours_at(v_max).size() == neighbours_at(neighbours.back()).size());

      std::size_t v_max = neighbours.back();
      //print(std::to_string(A[v].first) + " -> " + std::to_string(A[v_max].first), thread_id);
      std::size_t prev_max_clique_size = max_clique_size;
      branch_heuristic(v_max, vertex_neighbourhood(v_max, neighbours, max_clique_size),
                       max_clique_size, clique, depth + 1);
      if (max_clique_size > prev_max_clique_size) {
        // std::cout << clique.size() << "\n";
        clique.emplace_back(v);
      }
    }

    void branch_exact(std::size_t v, std::vector<std::size_t> &neighbours,
                      std::size_t &max_clique_size, std::vector<std::size_t> &clique,
                      std::size_t thread_id,
                      std::size_t depth = 1) const {
      if (neighbours.empty()) {
        if (depth > max_clique_size) {
          max_clique_size = depth;
          clique.clear();
          clique.reserve(depth);
          clique.emplace_back(v);
        }
        return;
      }

      while (not neighbours.empty() && depth + neighbours.size() > max_clique_size) {
        std::size_t u = neighbours.back();
        neighbours.pop_back();

        std::vector<std::size_t> new_neighbours = vertex_neighbourhood(u, neighbours, max_clique_size);
        std::size_t prev_max_clique_size = max_clique_size;

        auto begin = std::chrono::high_resolution_clock::now();

        branch_exact(u, new_neighbours, max_clique_size, clique, thread_id, depth + 1);

        auto end = std::chrono::high_resolution_clock::now();
        const double duration = std::chrono::duration<double>(end - begin).count();
        if (duration > .5 && depth == 1) {
          std::unique_lock print_lock(print_mtx);
          std::cerr << std::string((thread_id + 1), ' ') <<  thread_id << " branch of " << A[u].first << " took " << duration << "s\n";
        }

        if (prev_max_clique_size < max_clique_size) {
          clique.emplace_back(v);
        }
      }
    }

    inline std::vector<std::size_t>
    vertex_neighbourhood(std::size_t w, const std::vector<std::size_t> &neighbours,
                         const std::size_t max_clique_size) const {
      std::vector<std::size_t> new_neighbours;
      new_neighbours.reserve(neighbours_at(w).size());
      for (const auto &v : neighbours_at(w)) {
        if (neighbours_at(v).size() > max_clique_size) {
          for (const auto &u : neighbours) {
            if (v == u) {
              new_neighbours.emplace_back(u);
            }
          }
        }
      }

      // std::vector<std::size_t> tmp = new_neighbours;
      // for (const auto &v : new_neighbours) {
      //   //        std::cout << A[v].first << " " << A[v].second.size() << "\n";
      //   for (std::size_t u = 1; u < A[v].second.size(); ++u) {
      //     assert(A[u].second.size() <= A[u - 1].second.size());
      //   }
      // }

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

  void print_clique(const std::vector<vertex_t> &mc) {
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
    std::this_thread::sleep_for(3s);

    for (std::size_t num_turns = ac != 1 ? std::atol(av[1]) : -1; num_turns; num_turns--) {
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
