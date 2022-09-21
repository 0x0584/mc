#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <execution>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <memory_resource>
#include <numeric>
#include <set>
#include <thread>
#include <unordered_map>

#ifndef THREADS_PER_CORE
#define THREADS_PER_CORE 8
#endif

std::uint32_t num_threads =
    std::thread::hardware_concurrency() * THREADS_PER_CORE;

template <typename Func>
static inline auto benchmark(Func test_func, int iterations) {
  const auto start = std::chrono::system_clock::now();
  while (iterations > 0) {
    std::vector<std::thread> ths;
    ths.reserve(num_threads);
    std::uint32_t th_id = 0;
    while (iterations-- > 0 && th_id++ < num_threads) {
      ths.emplace_back(test_func);
    }
    for (auto &th : ths) {
      th.join();
    }
  }
  const auto stop = std::chrono::system_clock::now();
  const auto secs = std::chrono::duration<double>(stop - start);
  return secs.count();
}

template <typename Vec, typename Mtx>
static inline void colour(Vec &vec, Mtx &mtx, std::size_t num_v) {
  // vec.reserve(num_v);
  size_t num_col = 0;
  for (size_t i = 0; i < num_v; ++i) {
    std::size_t j = 0;
    while (std::any_of(vec[j].begin(), vec[j].end(),
                       [&](std::size_t k) { return mtx[i][k]; })) {
      j++;
    }
    vec[j].push_back(i);
    num_col = std::max(num_col, j + 1);
  }

  std::vector<size_t> cols;
  cols.reserve(num_v);
  std::vector<size_t> v;
  v.reserve(num_v);

  for (size_t col = 0; col < num_col; ++col) {
    std::fill_n(std::back_inserter(cols), vec[col].size(), col + 1);
    std::copy(vec[col].begin(), vec[col].end(), std::back_inserter(v));
  }
}

template <typename Vec, typename Mtx>
static inline void colour_par(Vec &vec, Mtx &mtx, std::size_t num_v) {
  // vec.reserve(num_v);
  size_t num_col = 0;
  for (size_t i = 0; i < num_v; ++i) {
    std::size_t j = 0;
    while (std::any_of(std::execution::par, vec[j].begin(), vec[j].end(),
                       [&](std::size_t k) { return mtx[i][k]; })) {
      j++;
    }
    vec[j].push_back(i);
    num_col = std::max(num_col, j + 1);
  }

  std::vector<size_t> cols;
  cols.reserve(num_v);
  std::vector<size_t> v;
  v.reserve(num_v);

  for (size_t col = 0; col < num_col; ++col) {
    std::fill_n(std::execution::par, std::back_inserter(cols), vec[col].size(),
                col + 1);
    std::copy(std::execution::par, vec[col].begin(), vec[col].end(),
              std::back_inserter(v));
  }
}

int main(int ac, const char *av[]) {
  int iterations{ac != 1 ? atoi(av[1]) : 500};

  std::cout << "iterations: " << iterations << " ";

  size_t num_v, num_e, num_exp;
  std::cin >> num_v >> num_e >> num_exp;

  std::unordered_map<size_t, std::set<size_t>> tmp_adj_lst;
  tmp_adj_lst.max_load_factor(0.5);
  tmp_adj_lst.reserve(num_v);
  for (size_t n_edges = num_e; n_edges--;) {
    size_t u, v;
    std::cin >> u >> v;
    tmp_adj_lst[u].emplace(v);
    tmp_adj_lst[v].emplace(u);
  }

  std::vector<std::vector<bool>> mtx(num_v, std::vector<bool>(num_v, false));
  for (const auto &[v, adj] : tmp_adj_lst) {
    for (const auto &u : adj) {
      mtx[v - 1][u - 1] = true;
    }
  }

  std::cout << "num_v " << num_v << " num_e " << num_e << "\n";

  double orig[4];

  {
    auto default_std_alloc = [&mtx, &num_v] {
      std::vector<std::vector<size_t>> vec(num_v);
      // for (auto &v : vec) {
      //   v.reserve(total_nodes);
      // }
      colour(vec, mtx, num_v);
    };

    auto default_pmr_alloc = [&mtx, &num_v] {
      std::pmr::vector<std::pmr::vector<size_t>> vec(num_v);
      // for (auto &v : vec) {
      //   v.reserve(total_nodes);
      // }
      colour(vec, mtx, num_v);
    };

    auto pmr_alloc_no_buf = [&mtx, &num_v] {
      std::pmr::monotonic_buffer_resource mbr;
      std::pmr::polymorphic_allocator<size_t> pa{&mbr};
      std::pmr::vector<std::pmr::vector<size_t>> vec(
          num_v, std::pmr::vector<size_t>(pa), pa);
      // for (auto &v : vec) {
      //   v.reserve(total_nodes);
      // }
      colour(vec, mtx, num_v);
    };

    auto pmr_alloc_and_buf = [&mtx, &num_v] {
      const auto buffer_size = num_v * sizeof(std::size_t);
      auto buffer = std::make_unique<std::byte[]>(buffer_size);
      std::pmr::monotonic_buffer_resource mbr{buffer.get(), buffer_size};
      std::pmr::polymorphic_allocator<int> pa{&mbr};
      std::pmr::vector<std::pmr::vector<int>> vec(
          num_v, std::pmr::vector<int>(pa), pa);
      // for (auto &v : vec) {
      //   v.reserve(total_nodes);
      // }
      colour(vec, mtx, num_v);
    };

    std::cout << std::fixed << std::setprecision(3);

    const double t1 = orig[0] = benchmark(default_std_alloc, iterations);
    std::cout << "t1 (default std alloc): " << t1 << " sec; t1/t1: " << t1 / t1
              << '\n';

    const double t2 = orig[1] = benchmark(default_pmr_alloc, iterations);
    std::cout << "t2 (default pmr alloc): " << t2 << " sec; t1/t2: " << t1 / t2
              << '\n';

    const double t3 = orig[2] = benchmark(pmr_alloc_no_buf, iterations);
    std::cout << "t3 (pmr alloc  no buf): " << t3 << " sec; t1/t3: " << t1 / t3
              << '\n';

    const double t4 = orig[3] = benchmark(pmr_alloc_and_buf, iterations);
    std::cout << "t4 (pmr alloc and buf): " << t4 << " sec; t1/t4: " << t1 / t4
              << '\n';
  }

  {
    auto default_std_alloc = [&mtx, &num_v] {
      std::vector<std::vector<size_t>> vec(num_v);
      // for (auto &v : vec) {
      //   v.reserve(total_nodes);
      // }
      colour_par(vec, mtx, num_v);
    };

    auto default_pmr_alloc = [&mtx, &num_v] {
      std::pmr::vector<std::pmr::vector<size_t>> vec(num_v);
      // for (auto &v : vec) {
      //   v.reserve(total_nodes);
      // }
      colour_par(vec, mtx, num_v);
    };

    auto pmr_alloc_no_buf = [&mtx, &num_v] {
      std::pmr::monotonic_buffer_resource mbr;
      std::pmr::polymorphic_allocator<size_t> pa{&mbr};
      std::pmr::vector<std::pmr::vector<size_t>> vec(
          num_v, std::pmr::vector<size_t>(pa), pa);
      // for (auto &v : vec) {
      //   v.reserve(total_nodes);
      // }
      colour_par(vec, mtx, num_v);
    };

    auto pmr_alloc_and_buf = [&mtx, &num_v] {
      const auto buffer_size = num_v * sizeof(std::size_t);
      auto buffer = std::make_unique<std::byte[]>(buffer_size);
      std::pmr::monotonic_buffer_resource mbr{buffer.get(), buffer_size};
      std::pmr::polymorphic_allocator<int> pa{&mbr};
      std::pmr::vector<std::pmr::vector<int>> vec(
          num_v, std::pmr::vector<int>(pa), pa);
      // for (auto &v : vec) {
      //   v.reserve(total_nodes);
      // }
      colour_par(vec, mtx, num_v);
    };

    std::cout << std::fixed << std::setprecision(3);

    const double t1 = benchmark(default_std_alloc, iterations);
    std::cout << "t1 (default std alloc): " << t1 << " sec; t1/t1: " << t1 / t1
              << " over " << orig[0] / t1 << '\n';

    const double t2 = benchmark(default_pmr_alloc, iterations);
    std::cout << "t2 (default pmr alloc): " << t2 << " sec; t1/t2: " << t1 / t2
              << " over " << orig[1] / t2 << '\n';

    const double t3 = benchmark(pmr_alloc_no_buf, iterations);
    std::cout << "t3 (pmr alloc  no buf): " << t3 << " sec; t1/t3: " << t1 / t3
              << " over " << orig[2] / t3 << '\n';

    const double t4 = benchmark(pmr_alloc_and_buf, iterations);
    std::cout << "t4 (pmr alloc and buf): " << t4 << " sec; t1/t4: " << t1 / t4
              << " over " << orig[3] / t4 << '\n';
  }
}
