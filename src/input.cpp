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

static inline char const *last_delimiter(char const *__restrict base,
                                         char const *__restrict tail) noexcept {
#define STRINGIFY(x) #x
#define PRAGMA_UNROLL(x) _Pragma(STRINGIFY(unroll x))
#define BLOCK_SCAN(BLOCK_SIZE)                                                 \
  do {                                                                         \
    while ((tail - base) >= BLOCK_SIZE) {                                      \
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

feed::buffer_content feed::read_chunk() {
  buffer_content chunk;
  char *base = chunk.buff.data();

  std::size_t total_sz = static_cast<std::size_t>(tail_rem - begin_rem);
  if (total_sz) {
    std::memcpy(base, begin_rem, total_sz);
  }

  const std::size_t to_read = CHUNK_SIZE - total_sz;
  in->read(base + total_sz, static_cast<std::streamsize>(to_read));
  if (in->bad()) [[unlikely]] {
    logger::error("FAILURE: cannot read from stream");
  }

  const std::size_t size_read = static_cast<std::size_t>(in->gcount());
  total_sz += size_read;

  char *tail = base + total_sz;
  char const *deli = last_delimiter(base, tail);
  if (!deli) {
    logger::error("FAILURE: chunk_size exceeded! recompile with a bigger size");
  }

  chunk.size = static_cast<std::size_t>(deli - base);

  const char *after_nl = (deli < tail) ? deli + 1 : tail;
  const std::size_t spill = static_cast<std::size_t>(tail - after_nl);
  if (spill) {
    std::memcpy(begin_rem, after_nl, spill);
  }
  tail_rem = begin_rem + spill;

  return chunk;
}
} // namespace mc
