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

feed::feed(input_source &in, std::size_t num_jobs)
    : in(in), batch_size(feed_jobs_scale(num_jobs)),
      remaining(memory::allocate<char>(BUFF_SIZE)), begin_rem(remaining),
      tail_rem(remaining), buffs(batch_size), read_idxs(batch_size),
      chunks(batch_size, {nullptr, nullptr, -1ul, this}) {
  for (auto &buff : buffs) {
    buff = memory::allocate<char>(BUFF_SIZE);
  }
  std::iota(read_idxs.begin(), read_idxs.end(), 0ul);
  fetch_idxs.reserve(batch_size);
  reader = std::jthread([this] mutable {
    std::pmr::vector<std::size_t> batch(memory::pool());
    batch.reserve(batch_size);
    while (true) {
      {
        std::unique_lock lk(recycle_mtx);
        recycle_cv.wait(lk,
                        [this] { return !read_idxs.empty() || !reading(); });
        if (!reading()) [[unlikely]] {
          std::size_t size = static_cast<std::size_t>(tail_rem - begin_rem);
          if (size != 0) {
            Assert(!read_idxs.empty());
            std::size_t idx = read_idxs.back();
            buffer buff = buffs[idx];

            std::memcpy(buff, begin_rem, size);

            chunks[idx].begin = buff;
            chunks[idx].end = buff + size;
            chunks[idx].idx = idx;

            tail_rem = begin_rem;

            {
              std::scoped_lock fk(fetch_mtx);
              fetch_idxs.push_back(idx);
              fetch_cv.notify_one();
            }
          }
          {
            std::scoped_lock fk(fetch_mtx);
            done.store(true, std::memory_order_release);
            fetch_cv.notify_all();
          }
          break;
        }
        logger::warn("read batch", logger::number_unit(read_idxs.size()));
        batch.insert(batch.end(), read_idxs.begin(), read_idxs.end());
        read_idxs.clear();
      }
#pragma unroll 8
      for (std::size_t i = 0; i < batch.size(); i++) {
        if (reading()) [[likely]] {
          read_chunk(batch[i]);
        } else {
          break;
        }
      }
      batch.clear();
    }
  });

  logger::debug("Feed has", feed_jobs_scale(num_jobs), "buffers for", num_jobs,
                "jobs");

  const std::size_t estimated_chunks = 1 + args::stream_size() / BUFF_SIZE;
  const double avg_edge_line =
      static_cast<double>(args::stream_size()) / in.num_e;
  const std::size_t edges_per_chunk = static_cast<std::size_t>(
      std::max(1.0, std::floor(BUFF_SIZE / avg_edge_line)));
  logger::info("Stream Size", logger::size_unit(args::stream_size()),
               "and Chunk Size", logger::size_unit(BUFF_SIZE));
  logger::info("Estimating", logger::number_unit(estimated_chunks),
               "Chunks with", logger::number_unit(edges_per_chunk),
               "Edges per Chunk");
}

feed::~feed() {
  memory::deallocate(remaining);
  for (auto &buff : buffs) {
    memory::deallocate(buff);
  }
  if (tail_rem != begin_rem) {
    // FIXME: use either exceptions or logger::error
    logger::error("INVALID file: no NL at the end of the file");
  }
  logger::debug("~feed()");
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

feed::chunk feed::fetch_chunk() {
  thread_local auto last = std::chrono::high_resolution_clock::now();

  std::size_t idx;
  {
    std::unique_lock lk(fetch_mtx);
    fetch_cv.wait(lk, [this] {
      return !fetch_idxs.empty() || done.load(std::memory_order_acquire);
    });
    if (fetch_idxs.empty() && done.load(std::memory_order_acquire))
        [[unlikely]] {
      return {nullptr, nullptr, 0, this};
    }
    idx = fetch_idxs.back();
    fetch_idxs.pop_back();
  }

  auto now = std::chrono::high_resolution_clock::now();
  chunks[idx].fetch_t = now - last;
  last = now;

  return chunks[idx];
}

void feed::read_chunk(std::size_t idx) {
  thread_local auto last = std::chrono::high_resolution_clock::now();

  Assert(tail_rem >= begin_rem);
  std::size_t total_sz = static_cast<std::size_t>(tail_rem - begin_rem);
  Assert(total_sz <= BUFF_SIZE);

  buffer base = buffs[idx];
  if (total_sz) {
    std::memcpy(base, begin_rem, total_sz);
  }

  const std::size_t to_read = BUFF_SIZE - total_sz;
  in->read(base + total_sz, static_cast<std::streamsize>(to_read));
  if (in->bad()) {
    logger::error("FAILURE: cannot read from stream");
  }

  total_sz += static_cast<std::size_t>(in->gcount());
  char *tail = base + total_sz;

  const char *deli = last_delimiter(base, tail);
  if (!deli) {
    throw std::runtime_error("Record exceeds BUFF_SIZE");
  }

  const std::size_t size = static_cast<std::size_t>(deli - base);
  const char *after_nl = (deli < tail) ? deli + 1 : tail;
  const std::size_t spill = static_cast<std::size_t>(tail - after_nl);
  if (spill) {
    std::memcpy(begin_rem, after_nl, spill);
  }
  tail_rem = begin_rem + spill;

  chunks[idx].begin = base;
  chunks[idx].end = base + size;
  chunks[idx].idx = idx;

  auto now = std::chrono::high_resolution_clock::now();
  chunks[idx].read_t = now - last;
  last = now;

  {
    std::scoped_lock lk(fetch_mtx);
    fetch_idxs.push_back(idx);
    fetch_cv.notify_one();
  }
}

} // namespace mc
