#include <iostream>
#include <set>
#include <unordered_map>
#include <vector>

int main() {
  std::size_t num_v, num_e, expected_max_clique_size;

  std::cin >> num_v >> num_e >> expected_max_clique_size;
  std::unordered_map<std::size_t, std::set<std::size_t>> tmp_adj_lst;

  tmp_adj_lst.max_load_factor(0.5);
  tmp_adj_lst.reserve(num_v);
  for (size_t n_edges = num_e; n_edges--;) {
    std::size_t u, v;
    std::cin >> u >> v;
    tmp_adj_lst[u].emplace(v);
    tmp_adj_lst[v].emplace(u);
  }

  std::vector<std::size_t> v_keys;
  v_keys.reserve(num_v);

  std::vector<std::vector<std::size_t>> A;
  std::vector<std::vector<bool>> M;
  M.resize(num_v);

  for (const auto &[v, N] : tmp_adj_lst) {
    v_keys.emplace_back(v);
  }
}
