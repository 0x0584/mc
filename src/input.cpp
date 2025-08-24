#include "input.hpp"
#include "core.hpp"

namespace mc {
std::uint16_t args::num_threads = thread::num_available_threads;

void args::parse(int argc, char *argv[]) {
  num_threads = thread::num_available_threads;
  logger::warn(thread::num_available_threads);
  for (int ch; (ch = getopt(argc, argv, "r:i:s:u:l:deyh")) != -1;) {
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

  logger::info("Running", args::exec_mode);

  if (stdin) {
    logger::info("Reading from STDIN");
  } else {
    logger::info("Reading from", filename);
  }
}

std::string feed::read_chunk() {
  buffer buff;
  const auto size_remaining = std::distance(remaining.begin(), tail_remaining);
  std::move(remaining.begin(), tail_remaining, buff.data());

  auto read = read_next(CHUNK_SIZE - size_remaining);
  auto begin_read = read.first.begin();
  const auto size_read = read.second;
  std::move(begin_read, begin_read + size_read, buff.data() + size_remaining);

  const auto tail_buff = buff.begin() + (size_remaining + size_read);
  auto delimiter = tail_buff;
  while (delimiter != buff.begin() && *--delimiter != feed::deli)
    ;
  if (delimiter == buff.begin()) {
    logger::error("FAILURE: chunk_size=", CHUNK_SIZE,
                  " exceeded! recompile with a bigger size");
  }

  const auto buffer_size =
      static_cast<std::size_t>(std::distance(buff.begin(), delimiter++));
  std::move(delimiter, tail_buff, remaining.data());
  tail_remaining = remaining.begin() + std::distance(delimiter, tail_buff);

  return std::string{std::move(buff.data()), buffer_size};
}
} // namespace mc
