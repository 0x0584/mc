#include "include/graph.hpp"

struct max_clique {
  using vertex_t = unsigned;
  using edge_t = std::pair<vertex_t, vertex_t>;

  mutable std::unordered_map<vertex_t, std::unordered_map<vertex_t, bool>> edges;

  max_clique(const std::unordered_set<edge_t> &E) : A(adjacency_list(E)), D(degree_count()),
                                                    Kore(core_numbers()), H(heuristic()) {
    prune_graph();
    std::cout << "Current Clique of size " << H.size() << ": ";
        for (const auto &u : H) {
          std::cout << u << " ";
        }
        std::cout << "\n";

    std::size_t size = H.size();
    while (not D.empty()) {
      if (size < H.size()) {
        std::cout << "Current Clique: " << H.size() << ": ";
        for (const auto &u : H) {
          std::cout << u << " ";
        }
        std::cout << "\n";
        size = H.size();
      }
      //std::cout << "Graph has: ";
      std::set<vertex_t> S;
      for (const auto &v : A) {
        S.emplace(v.first);
      }
      // for (const auto &v : S) {
      //   //std::cout << v << "(" << A.at(v).size() << ") ";
      // }
      //std::cout << "\n";

      auto it = D.begin();
      vertex_t v = D.begin()->second;
      it++;
      //std::cout << "Picking vertex: " << v << " with degree " << D.begin()->first << "\n";
      // if (it != D.end()) {
      //   //std::cout << "Picking vertex: " << it->second << " with degree " << it->first << "\n";
      // }
      std::unordered_set<vertex_t> P = neighbourhood(v);

      //std::cout << " Has: ";
      // for (const auto &e : P) {
      //   //std::cout << e << " ";
      // }
      //std::cout << "\n";

      prune_neighbourhood(P);

      //std::cout << " After prune Has: ";
      // for (const auto &e : P) {
      //   //std::cout << e << " ";
      // }
      //std::cout << "\n";

      if (num_greedy_colours(P) > H.size()) {
        branch(std::move(P));
      }
      prune_vertex(v);
      //std::cout << " --------------------- Looping ----------------------\n";
    }
  }

  void prune_neighbourhood(std::unordered_set<vertex_t> &P) const {
    if (P.size() <= H.size()) {
      //std::cout << " But " << P.size() << " < " << H.size() << " clique size\n";
      P = {};
      return;
    }

    std::multimap<std::size_t, vertex_t> K;
    core_numbers(K, P);
    std::size_t K_max = K.rbegin()->first + 1;
    if (K_max <= H.size()) {
      //std::cout << " But Largest degree " << K_max << " < " << H.size() << " clique size\n";
      P = {};
      return;
    }

    for (const auto &v : K) {
      if (v.first <= H.size()) {
        //std::cout << "  Pruning " << v.second << "\n";
        P.erase(v.second);
      }
    }
  }

  inline std::unordered_set<vertex_t> neighbourhood(vertex_t v) const {
    std::unordered_set<vertex_t> P;
    P.emplace(v);
    for (const auto &u : A.at(v)) {
      if (Kore.at(u) > H.size()) {
          P.emplace(u);
        if (!edges[u][v]) {

          //        assert(edges[u][v] == false);
          edges[u][v] = true;
        }
      }
    }
    return std::move(P);
  }

  inline std::unordered_set<vertex_t>
  induce_neighbourhood(const std::unordered_set<vertex_t> &P, vertex_t u) const {
    std::unordered_set<vertex_t> Q = neighbourhood(u);
    //std::cout << " neighbourhood of " << u << ": ";
    // for (const auto &e : Q) {
    //   //std::cout << e << " ";
    // }
    //std::cout << "\n";
    // std::vector<vertex_t> R;
    // std::set_intersection(P.begin(), P.end(), Q.begin(), Q.end(), std::back_inserter(R));

    std::unordered_set<vertex_t> R;

    const auto &small = P.size() > Q.size() ? Q : P;
    const auto &big = P.size() > Q.size() ? P : Q;
    for (const auto &v : small) {
      if (big.count(v)) {
        R.emplace(v);
      }
    }
    return std::move(R);
  }

  static inline std::unordered_set<vertex_t>
  enlarge_clique(const std::unordered_set<vertex_t> &C, vertex_t u) {
    std::unordered_set<vertex_t> C2 = C;
    C2.emplace(u);
    return std::move(C2);
  }

  void branch(std::unordered_set<vertex_t> P, std::unordered_set<vertex_t> C = {}, int depth = 1) {
    // auto print =
    //   [&depth]() -> std::ostream & {
    //     return std::cout << std::string(depth * 1.5, ' ');
    //   };

    // print() << " > Branch at depth " << depth << " with ";
    // for (const auto &e : P) {
    //   std::cout << e << " ";
    // }
    // print() << "\n";
    // if (not C.empty()) {
    //   print() << "  > Best clique so far ";
    //   for (const auto &e : C) {
    //     std::cout << e << " ";
    //   }
    //   print() << "\n";
    // }

    while (not P.empty() && P.size() + C.size() > H.size()) {
      vertex_t u = *P.begin();
      P.erase(u);

      //   print() << "   Picking " << u << "\n";

      std::unordered_set<vertex_t> C2 = enlarge_clique(C, u);
      std::unordered_set<vertex_t> P2 = induce_neighbourhood(P, u);

      // print() << "   induced nei ";
      // for (const auto &e : P2) {
      //   std::cout << e << " ";
      // }
      // print() << "\n";

      // print() << "   Enlarged Clique ";
      // for (const auto &e : C2) {
      //   //std::cout << e << " ";
      // }
      // print() << "\n";

      if (not P2.empty()) {
        if (C2.size() + num_greedy_colours(P2) > H.size()) {
          branch(std::move(P2), std::move(C2), depth + 1);
        } else {
          //          print() << " ### Too small clique\n";
        }
      } else if (C2.size() > H.size()) {
        // std::cout << "Found better Clique\n";
        // std::cout << "  > Best clique so far ";
        // for (const auto &e : C) {
        //   std::cout << e << " ";
        // }
        // std::cout << "\n";
        // char c;
        // std::cin >> c;
        H = std::move(C2);
        prune_graph();
      } else {
        //        print() << " ### clique not found \n";
      }


    }
  }

  std::size_t num_greedy_colours(const std::unordered_set<vertex_t> &S) const {
    if (S.empty()) {
      return 0;
    }

    std::set<std::size_t> C;    // colours
    std::unordered_map<vertex_t, std::size_t> L; // coloured vertices
    std::multimap<std::size_t, vertex_t> K;
    core_numbers(K, S);
    for (const auto &v : K) {
      const auto &A_v = A.at(v.second);
      std::set<std::size_t> N_C; // neighbour colours
      std::vector<vertex_t> N;

      for (const auto &u : A_v) {
        if (S.count(u)) {
          //  //std::cout << " Vertex " << u << " has colour " << L[u] << "\n";
          N_C.emplace(L[u]);
        }
      }
      N_C.erase(0);             // fresh vertices have colour 0

      std::size_t colour = 1;
      for (const auto &c : N_C) {
        if (c != colour) {
          break;
        }
        colour++;
      }
      L[v.second] = colour;
      //      //std::cout << " Vertex : " << v.second << " has colour " << colour << "\n";
      C.emplace(colour);
    }

    // //std::cout << "Colours\n";
    // for (const auto &l : L) {
    //   //std::cout << " vertex " << l.first << " has colour " << l.second << "\n";
    // }

    // //std::cout << " >> Colours: ";
    // for (const auto &c : C) {
    //   //std::cout << c << "\n";
    // }
    //std::cout << " Num Colours " << C.size() << "\n";
    return C.size();
  }

  std::unordered_map<vertex_t, std::unordered_set<vertex_t>>
  adjacency_list(const std::unordered_set<edge_t> &E) {
    std::unordered_map<vertex_t, std::unordered_set<vertex_t>> adj;
    //std::cout << "Adjacency List\n";
    for (const auto &e : E) {
      adj[e.first].emplace(e.second);
      adj[e.second].emplace(e.first);
    }

    edges.reserve(adj.size());
    for (const auto &v : adj) {
      edges[v.first].reserve(v.second.size());
      //std::cout << "Vertex: " << v.first << " { ";
      std::unordered_map<vertex_t, bool> edge;
      for (const auto &u : v.second) {
        edges[v.first][u] = false;
      }

      //      edges.emplace(v, edge);
      //std::cout << "}\n";
    }

    return std::move(adj);
  }

  std::multimap<std::size_t, vertex_t> degree_count() {
    //std::cout << "Degree Count\n";
    std::multimap<std::size_t, vertex_t> degree_lst;
    for (const auto &v : A) {
      degree_lst.emplace(v.second.size(), v.first);
    }
    return std::move(degree_lst);
  }

  std::unordered_set<vertex_t> heuristic() {
    //std::cout << "Heuristic Count\n";
    std::multimap<std::size_t, vertex_t,
                  std::greater<std::size_t>> K; // core numbers of all vertices
    core_numbers(K);

    std::unordered_set<vertex_t> H;
    for (const auto &v : K) {
      if (v.first < H.size()) {
        break;
      }

      //      //std::cout << "vertex " << v.second << " has K " << v.first << "\n";

      std::multimap<std::size_t, vertex_t,
                    std::greater<std::size_t>> S; // core number of v neighbourhood

      core_numbers(S, A.at(v.second));

      // if (not S.empty()) {
      //   //std::cout << "\n > S neighbours of " << v.second << " { ";
      //   for (const auto &e : S) {
      //     //std::cout << e.second << " ";
      //   }
      //   //std::cout << "}\n";
      // } else {
      //   //std::cout << " NO NEIGHBOPURS FOR "  << v.second << "\n";
      //      }

      std::unordered_set<vertex_t> C; // tmp clique
      C.emplace(v.second);

      for (const auto &u : S) {
        if (is_clique(C, u.second)) {
          C.emplace(u.second);
        }
      }

      // if (not C.empty()) {
      //   //std::cout << "\n > Tmp clique of " << v.second << " { ";
      //   for (const auto &e : C) {
      //     //std::cout << e << " ";
      //   }
      //   //std::cout << "}\n";
      // }

      if (C.size() > H.size()) {
        H = std::move(C);
      }
      //      //std::cout << "\n";
    }

    // //std::cout << "\n > Heuristic clique SIZE = " << H.size() << " { " ;
    // for (const auto &e : H) {
    //   //std::cout << e << " ";
    // }
    // //std::cout << "}\n";

    return std::move(H);
    // return {};
  }

  std::unordered_map<vertex_t, std::size_t> core_numbers() {
    //std::cout << "Core Numbers Count\n";
    std::unordered_map<std::size_t, vertex_t> W;
    W.reserve(A.size());
    std::unordered_map<std::size_t, std::unordered_set<std::size_t>> A2;
    A2.reserve(A.size());
    {
      std::unordered_map<vertex_t, std::size_t> W2;
      W2.reserve(A.size());
      std::size_t index = 1;
      const auto vertex_index =
        [&W, &W2, &index](const vertex_t &u) {
          auto &u_index = W2[u];
          if (u_index == 0) {
            W[u_index = index++] = u;
          }
          return u_index;
        };
      for (const auto &u : A) {
        auto &S_u = A2[vertex_index(u.first)];
        //      //std::cout << "vertex " << e.first << "=" << u << "\n";
        for (const auto &v : u.second) {
          S_u.emplace(vertex_index(v));
        }
        // auto v = vertex_index(e.second);
        // ////std::cout << "vertex " << e.second << "=" << v << "\n";
        // A[u].emplace(v);
        // A[v].emplace(u);
      }
    }

    std::vector<std::size_t> Deg(A.size() + 1);
    for (const auto &v : A2) {
      Deg[v.first] = v.second.size();
      ////std::cout << "Vertex " << W[v.first] << " Has Degree " << v.second.size() << "\n";
    }
    {
      std::size_t todo = A.size(), level = 1;
      const auto process_sublevel =
        [&Deg, &W, &A2, &todo](std::unordered_set<std::size_t> curr, std::size_t level) {
          std::unordered_set<std::size_t> next;
          do {
            todo -= curr.size();
            for (const auto &v : curr) {
              //            //std::cout << " > Vertex: " << W[v] << "\n";
              for (const auto &u : A2.at(v)) {
                //              //std::cout << "  --> Vertex: " << W[u];
                if (Deg[u] > level) {
                  if (--Deg[u] == level) {
                    //      //std::cout << " picked";
                    next.emplace(u);
                  }
                  //K[W[u]] = level;
                }
                //      //std::cout << "\n";
              }
            }
            ////std::cout << " ----------- \n";
            curr = std::move(next);
          } while (not curr.empty());
        };

      for (; todo > 0; ++level) {
        std::unordered_set<std::size_t> curr;
        // //std::cout << "scan level " << level << "\n";
        for (std::size_t i = 1; i < Deg.size(); ++i) {
          if (Deg[i] == level) {
            ////std::cout << "Picked " << W[i] << " with degree " << D[i] << "\n";
            curr.emplace(i);
            //K[W[i]] = level;
          }
        }
        //      //std::cout << "\n";

        if (not curr.empty()) {
          process_sublevel(std::move(curr), level);
        }
      }
    }

    std::unordered_map<vertex_t, std::size_t> K;
    for (std::size_t i = 1; i < Deg.size(); ++i) {
      K[W[i]] = Deg[i];
    }

    // for (std::size_t i = 1; i < D.size(); ++i) {
    //   //std::cout << "degree of " << W[i] << " is " << D[i];
    //   K[W[i]] = D[i];
    //   //std::cout << "\n";
    // }
    // //std::cout << "\n";


    // for (std::size_t i = 1; i < Deg.size(); ++i) {
      //std::cout << " > K[" << W[i] << "]=" << Deg[i] << "\n";
    // }

    return std::move(K);
  }

  template <typename Container>
  inline void core_numbers(Container &K) const {
    for (const auto &v : A) {
      K.emplace(Kore.at(v.first), v.first);
    }
  }

  template <typename DstContainer, typename SrcContainer>
  inline void core_numbers(DstContainer &K, const SrcContainer &S) const {
    for (const auto &u : S) {
      K.emplace(Kore.at(u), u);
    }
  }

  template <typename Callable>
  inline bool reduce_vertex_degree(vertex_t u, const std::size_t degree, Callable func) {
    //    //std::cout << "  Checking: " << u << " == \n";
    const auto &flag = func();
    for (auto u_itr = D.find(degree); u_itr != D.end(); u_itr++) {
      if (u_itr->second == u) {
        //  //std::cout << "  Found: " << u << " " << v->first << " == \n";
        D.erase(u_itr);
        if (flag) {
          D.insert({u_itr->first - 1, u});
        } else {
          assert(A.count(u));
          A.erase(u);
          assert(Kore.count(u));
          Kore.erase(u);
        }
        return flag;
      } else {
        //        //std::cout << v->first << " " << A.at(u).size() << "\n";
        assert(u_itr->first == degree);
      }
    }
    return false;
  }

  inline void prune_vertex(vertex_t v) {
    if (A.count(v)) {
      //std::cout << " >>>>>> PRUNE " << v << "\n";

      auto &S_v = A.at(v);
      for (const auto &u : S_v) {
        assert(A.count(u));
        auto &S_u = A.at(u);
        assert(not S_u.empty());
        const std::size_t deg = S_u.size();
        if (reduce_vertex_degree(u, deg, [&deg] { return deg > 1; })) {
          S_u.erase(v);
        } else {
          //std::cout << " |>>>>>> PRUNE " << u << "\n";
        }
      }
      reduce_vertex_degree(v, S_v.size(), [] { return false; });
    }
  }

  inline void prune_graph() {
    //std::cout << " >>>>>>>>> Cleaning Graph\n";
    std::unordered_set<vertex_t> K;
    for (const auto &v : Kore) {
      if (v.second < H.size()) {
        //std::cout << "Vertex " << v.first << " has core " << v.second
        //                  << " while Heuristic is " << H.size() << "\n";
        K.emplace(v.first);
      }
    }
    while (not K.empty()) {
      prune_vertex(*K.begin());
      K.erase(*K.begin());
    }
    Kore = std::move(core_numbers());
  }

  inline bool is_clique(const std::unordered_set<vertex_t> &C, vertex_t u) const {
    for (const auto &w : C) {
      const auto &N = A.at(w); // neighbors of w
      if (std::find(N.begin(), N.end(), u) == N.end()) {
        return false;
      }
    }
    return true;
  }

  const std::unordered_set<vertex_t> &operator()() const {
    return H;
  }

  std::unordered_set<vertex_t> operator()() {
    return std::move(H);
  }

  std::unordered_set<vertex_t>::const_iterator begin() const {
    return H.begin();
  }

  std::unordered_set<vertex_t>::const_iterator end() const {
    return H.end();
  }

private:
  std::unordered_map<vertex_t, std::unordered_set<vertex_t>> A; // Adjacency List
  std::multimap<std::size_t, vertex_t> D;         // Degree List
  std::unordered_map<vertex_t, std::size_t> Kore; // Core Numbers
  std::unordered_set<vertex_t> H;                 // heuristic clique
};

int main() {
  std::size_t num_v, num_e;
  std::cin >> num_v >> num_e;

  std::unordered_set<max_clique::edge_t> E;


  // E.emplace(1, 2);
  // E.emplace(1, 4);
  // E.emplace(1, 3);
  // E.emplace(1, 8);

  // E.emplace(2, 6);
  // E.emplace(2, 5);
  // //  E.emplace(2, 4);
  // E.emplace(2, 7);

  // E.emplace(3, 4);
  // E.emplace(3, 4);
  // E.emplace(3, 6);
  // E.emplace(3, 5);
  // E.emplace(3, 10);

  // E.emplace(4, 6);
  // E.emplace(4, 7);
  // E.emplace(4, 5);
  // E.emplace(4, 9);
  // E.emplace(4, 10);

  // E.emplace(6, 8);

  // E.emplace(5, 6);
  // E.emplace(5, 7);
  // E.emplace(5, 8);

  // E.emplace(7, 9);
  // E.emplace(7, 9);
  // E.emplace(7, 10);

  // E.emplace(8, 9);

  // E.emplace(9, 10);


  // for (const auto &e : E) {
    //std::cout << e.first << "-" << e.second << "\n";
  //  }

  //  std::vector<std::vector<bool>> mtx(num_v + 1, std::vector<bool>(num_v + 1));
  while (num_e--) {
    std::size_t u, v;
    std::cin >> u >> v;
    //mtx[u][v] = true;
    E.emplace(u, v);
  }

  // for (const auto &row : mtx) {
  //   for (const auto &col : row) {
  //     //std::cout << col << ", ";
  //   }
  //   //std::cout << "\n";
  // }
  // //std::cout << "\n";
  // return 0;
  max_clique mc(E);
  std::set<max_clique::vertex_t> m(mc.begin(), mc.end());
  std::cout << "SIZE = " << m.size() << " { ";
  for (const auto &e : m) {
    std::cout << e << " ";
  }
  std::cout << "}\n";
  return 0;
}
