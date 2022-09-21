#include <algorithm>
#include <assert.h>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using vertex_t = std::size_t;
using edge_t = std::pair<vertex_t, vertex_t>;

template <typename T>
inline void hash_combine(std::size_t &seed, T const &key) {
  std::hash<T> hasher;
  seed ^= hasher(key) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
template <typename T1, typename T2> struct hash<std::pair<T1, T2>> {
  std::size_t operator()(std::pair<T1, T2> const &p) const {
    std::size_t seed1(0);
    ::hash_combine(seed1, p.first);
    ::hash_combine(seed1, p.second);

    std::size_t seed2(0);
    ::hash_combine(seed2, p.second);
    ::hash_combine(seed2, p.first);

    return std::min(seed1, seed2);
  }
};
} // namespace std

std::map<vertex_t, std::size_t>
core_numbers(const std::unordered_set<edge_t> &E) {
  std::unordered_map<std::size_t, vertex_t> W;
  std::unordered_map<std::size_t, std::unordered_set<std::size_t>> A;
  {
    std::unordered_map<vertex_t, std::size_t> W2;
    std::size_t index = 1;
    const auto vertex_index = [&W, &W2, &index](const vertex_t &u) {
      auto &u_index = W2[u];
      if (u_index == 0) {
        W[u_index = index++] = u;
      }
      return u_index;
    };
    for (const auto &e : E) {
      auto u = vertex_index(e.first);
      //      std::cout << "vertex " << e.first << "=" << u << "\n";
      auto v = vertex_index(e.second);
      // std::cout << "vertex " << e.second << "=" << v << "\n";
      A[u].emplace(v);
      A[v].emplace(u);
    }
  }

  std::vector<std::size_t> D(A.size() + 1);
  for (const auto &v : A) {
    D[v.first] = v.second.size();
    // std::cout << "Vertex " << W[v.first] << " Has Degree " << v.second.size()
    // << "\n";
  }
  {
    std::size_t todo = A.size(), level = 1;
    const auto process_sublevel = [&D, &W, &A,
                                   &todo](std::unordered_set<std::size_t> curr,
                                          std::size_t level) {
      std::unordered_set<std::size_t> next;
      do {
        todo -= curr.size();
        for (const auto &v : curr) {
          //            std::cout << " > Vertex: " << W[v] << "\n";
          for (const auto &u : A.at(v)) {
            //              std::cout << "  --> Vertex: " << W[u];
            if (D.at(u) > level) {
              if (--D.at(u) == level) {
                //      std::cout << " picked";
                next.emplace(u);
              }
              // K[W[u]] = level;
            }
            //      std::cout << "\n";
          }
        }
        // std::cout << " ----------- \n";
        curr = std::move(next);
      } while (not curr.empty());
    };

    for (; todo > 0; ++level) {
      std::unordered_set<std::size_t> curr;
      // std::cout << "scan level " << level << "\n";
      for (std::size_t i = 1; i < D.size(); ++i) {
        if (D[i] == level) {
          // std::cout << "Picked " << W[i] << " with degree " << D[i] << "\n";
          curr.emplace(i);
          // K[W[i]] = level;
        }
      }
      //      std::cout << "\n";

      if (not curr.empty()) {
        process_sublevel(std::move(curr), level);
      }
    }
  }

  std::map<vertex_t, std::size_t> K;
  for (std::size_t i = 1; i < D.size(); ++i) {
    K[W[i]] = D[i];
  }

  // for (std::size_t i = 1; i < D.size(); ++i) {
  //   std::cout << "degree of " << W[i] << " is " << D[i];
  //   K[W[i]] = D[i];
  //   std::cout << "\n";
  // }
  // std::cout << "\n";

  // for (std::size_t i = 1; i < D.size(); ++i) {
  //   std::cout << " > K[" << W[i] << "]=" << D[i] << "\n";
  // }
  return std::move(K);
}

int main() {
  std::size_t num_v, num_e;
  std::cin >> num_v >> num_e;

  //  std::cout << "V=" << num_v << " E=" << num_e << "\n";
  std::unordered_set<edge_t> E;
  while (num_e--) {
    std::size_t u, v;
    std::cin >> u >> v;
    //    std::cout << u << "-" << v << "\n";
    E.emplace(u, v);
  }

  //  std::cout << E.size() << "\n";
  for (const auto &v : core_numbers(E)) {
    std::cout << "K(" << v.first << ")=" << v.second << "\n";
  }
  std::cout << "\n";
  return 0;
}
