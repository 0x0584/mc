#include "input.hpp"
#include "core.hpp"

#include <cstring>

namespace mc {
std::uint16_t args::num_threads = thread::num_available_threads;

void args::parse(int argc, char *argv[]) {
  num_threads = thread::num_available_threads;
  for (int ch; (ch = getopt(argc, argv, "t:r:i:s:u:l:deyh")) != -1;) {
    switch (ch) {
    case 'r':
      num_turns = std::max(1l, std::atol(optarg));
      break;
    case 'i':
      stdin = false;
      file = std::ifstream(filename = optarg);
      break;
    case 'd':
      undirected = false;
      break;
    case 's':
      expect_size = true;
      size = static_cast<std::size_t>(std::atol(optarg));
      logger::info("Expecting a Max Clique of size", size);
      break;
    case 'u':
      upper_bound = static_cast<std::size_t>(std::atol(optarg));
      logger::info("Expected Upper Bound for Max Clique of size", upper_bound);
      break;
    case 'l':
      lower_bound = static_cast<std::size_t>(std::atol(optarg));
      logger::info("Expected Lower Bound for Max Clique of size", lower_bound);
      break;
    case 'e':
      if (exec_mode == flavour::heuristic) {
        // just in case both hybrid and excat were specified, run as hybrid
        exec_mode = flavour::exact;
      }
      break;
    case 'y':
      exec_mode = flavour::hybrid;
      break;
    case 'o':
      draw = false;
      break;
    case 't': {
      int arg_num_threads = std::atoi(optarg);
      num_threads = static_cast<std::uint16_t>(arg_num_threads);
      if (arg_num_threads <= 0 || num_threads > thread::num_available_threads) {
        throw std::runtime_error("too many threads! abort.");
      }
      break;
    }
    case 'h':
    case ':':
    case '?':
    default:
      std::cerr << "Find the Max Clique of a graph using branch-and-bound\n\n"
                << "  -r N run for N times\n"
                << "  -t N number of threads\n"
                << "  -i FILE take input from FILE instead of STDIN\n"
                << "  -u N expect a at most a clique of size N\n"
                << "  -l N expect at least a clique of size N\n"
                << "  -e run the algorithm as EXACT (default HEURISTIC)\n"
                << "  -y run the algorithm as HYBRID (HEURISTIC + EXACT)\n"
                << "  -o output a Graphiz Dot file of the graph and the clique"
                << "  -d use DIRECTED edges instead of the default UNDIRECTED\n"
                << "\n";
      exit(EXIT_FAILURE);
    }
  }

  logger::info("Running", exec_mode, "with", num_threads, "Threads");

  if (stdin) {
    logger::info("Reading from STDIN");
  } else {
    logger::info("Reading from", filename);
  }
}

feed::feed(input_source &in)
    : in(in), _estimated_chunks(1 + args::stream_size() / BUFF_SIZE),
      remaining(BUFF_SIZE, memory::pool()), begin_rem(remaining.data()),
      tail_rem(remaining.data()),
      buffs(_estimated_chunks, buffer(memory::pool()), memory::pool()),
      was_read(_estimated_chunks) {
  const double avg_edge_line =
      static_cast<double>(args::stream_size()) / in.num_e;
  const std::size_t _edges_per_chunk = static_cast<std::size_t>(
      std::max(1.0, std::floor(BUFF_SIZE / avg_edge_line)));
  logger::info("Stream Size", logger::size_unit(args::stream_size()),
               "and Chunk Size", logger::size_unit(BUFF_SIZE));
  logger::info("Estimating", logger::number_unit(_estimated_chunks),
               "Chunks with", logger::number_unit(_edges_per_chunk),
               "Edges per Chunk");
}

feed::~feed() {
  if (tail_rem != begin_rem) {
    // FIXME: use either exceptions or logger::error
    logger::error("INVALID file: no NL at the end of the file");
  }
  logger::debug("~feed()");
}

void feed::prepare() {
  buffer_recycle = std::jthread([this] {
    std::pmr::vector<std::size_t> pending(buffs.size(), memory::pool());
    std::iota(pending.begin(), pending.end(), 0ul);
    while (!pending.empty()) {
      {
        std::unique_lock lk(mtx);
        recycle_cv.wait_for(lk, std::chrono::milliseconds(100));
      }
      std::size_t alive = pending.size();
      std::size_t old_alive = alive;
#pragma unroll 32
      for (std::size_t i = 0; i < alive;) {
        const std::size_t idx = pending[i];
        if (was_read[idx].exchange(false, std::memory_order_relaxed)) {
          buffs[idx] = buffer{memory::pool()};
          std::swap(pending[i], pending[--alive]);
          pending.pop_back();
        } else {
          ++i;
        }
      }
      if (old_alive > alive) [[unlikely]] {
        logger::debug("reclaimed:", (old_alive - alive), "buffer");
      }
    }
  });
}

static inline char const *last_delimiter(char const *__restrict base,
                                         char const *__restrict tail) noexcept {
#define STRINGIFY(x) #x
#define PRAGMA_UNROLL(x) _Pragma(STRINGIFY(unroll x))
#define BLOCK_SCAN(BLOCK_SIZE)                                                 \
  do {                                                                         \
    if ((tail - base) >= BLOCK_SIZE) {                                         \
      tail -= BLOCK_SIZE;                                                      \
      PRAGMA_UNROLL(BLOCK_SIZE)                                                \
      for (int i = BLOCK_SIZE - 1; i >= 0; --i) {                              \
        if (tail[i] == feed::deli)                                             \
          return tail + i;                                                     \
      }                                                                        \
    }                                                                          \
  } while (0)

  BLOCK_SCAN(32);
  BLOCK_SCAN(16);
  BLOCK_SCAN(8);

  while (tail != base) {
    --tail;
    if (*tail == feed::deli)
      return tail;
  }

  return nullptr;
#undef STRINGIFY
#undef PRAGMA_UNROLL
#undef BLOCK_SCAN
}

feed::chunk feed::read_chunk() {
  if (!reading()) [[unlikely]] {
    logger::error("failure to read");
  }

  chunk fresh(buffs[current_buff], was_read[current_buff]);
  current_buff++;

  if (current_buff > buffs.size()) [[unlikely]] {
    logger::error("chunk calculation is not correct");
  }

  fresh.resize();

  char *base = fresh.begin();

  std::size_t total_sz = static_cast<std::size_t>(tail_rem - begin_rem);
  if (total_sz) [[likely]] {
    std::memcpy(base, begin_rem, total_sz);
  }

  const std::size_t to_read = BUFF_SIZE - total_sz;
  in->read(base + total_sz, static_cast<std::streamsize>(to_read));
  if (in->bad()) [[unlikely]] {
    logger::error("FAILURE: cannot read from stream");
  }

  const std::size_t size_read = static_cast<std::size_t>(in->gcount());
  total_sz += size_read;

  char *tail = base + total_sz;
  char const *deli = last_delimiter(base, tail);
  if (!deli) [[unlikely]] {
    logger::error("FAILURE:", BUFF_SIZE,
                  "exceeded! recompile with a bigger size");
  }

  fresh.resize(static_cast<std::size_t>(deli - base));

  const char *after_nl = (deli < tail) ? deli + 1 : tail;
  const std::size_t spill = static_cast<std::size_t>(tail - after_nl);
  if (spill) [[likely]] {
    std::memcpy(begin_rem, after_nl, spill);
  }
  tail_rem = begin_rem + spill;

  return fresh;
}
} // namespace mc
