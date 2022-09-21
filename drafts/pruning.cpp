#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>

std::size_t num_colours_greedy(const std::vector<std::size_t> &K,
                               std::vector<std::size_t> &S,
                               std::size_t max_clique_size) const {
  std::unordered_set<std::size_t> C;              // colours
  std::unordered_map<std::size_t, std::size_t> L; // vertex colour

  L.reserve(S.size());
  C.reserve(S.size());

  std::sort(
      S.begin(), S.end(),
      [&](const std::size_t &u, const std::size_t &v) { return K[u] < K[v]; });

  // for (const auto &v : S) {
  //   print(std::to_string(v));
  // }

  for (const auto &v : S) {
    // std::size_t v = *it;
    // print(std::string(" colour of ") + std::to_string(v));

    std::set<std::size_t> N_C; // neighbour colours
    for (const auto &v : neighbours_at(v)) {
      if (K[v] >= max_clique_size) {
        for (const auto &u : S) {
          if (u == v) {
            N_C.emplace(L[u]);
            break;
          }
        }
      }
    }
    N_C.erase(0); // fresh vertices have colour 0

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

  assert(C.size() <= S.size());
  return C.size();
}

std::vector<std::size_t> core_numbers(std::vector<std::size_t> K) {
  for (std::size_t todo = K.size(), level = 1; todo; ++level) {
    std::vector<std::size_t> curr;
    curr.reserve(K.size());

    for (std::size_t v = 0; v < K.size(); ++v) {
      if (K[v] == level) {
        curr.emplace_back(v);
      }
    }

    for (std::vector<std::size_t> next; not curr.empty();
         curr = std::move(next)) {
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

  for (auto &[_, neighbours] : A) {
    std::sort(neighbours.begin(), neighbours.end(),
              [&](const std::size_t &u, const std::size_t &v) {
                return K[u] < K[v];
              });
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
