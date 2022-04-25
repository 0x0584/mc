#ifndef MAX_CLIQUE_HPP
#define MAX_CLIQUE_HPP

#include "algorithm.hpp"
#include "reader.hpp"

namespace mc {
class max_clique {
  const reader &reader;

public:
  const size_type expected_size;

  max_clique() = delete;
  max_clique(const max_clique &) = delete;

  max_clique(const struct reader &reader, size_type expected_size = 0)
      : reader(reader), expected_size(expected_size) {
#ifndef NDEBUG
    std::cerr << reader;
#endif
  }

  max_clique &operator=(const max_clique &) = delete;

  using algorithm_type = algorithm_exact;
  using bound_strategy_type = bound_by_colouring;
  reader::vertex_set solve(std::uint32_t nthreads = num_threads) const {
    source(reader.get());
    source_clique source;

    std::shared_mutex thread_mtx;
    std::condition_variable_any pool;
    std::vector<std::thread> threads(nthreads);
    std::vector<bool> available(nthreads, true);

    bound_strategy_type branch_bound();
    while () {
    }

    reader::vertex_set clique;

    return clique;
  }
};
} // namespace mc

#endif // MAX_CLIQUE_HPP
