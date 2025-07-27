// input.hpp
//
// Copyright (C) 2024  0x0584 (Anas)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
// USA.

#ifndef INPUT_HPP
#define INPUT_HPP

#include <cstdlib>
#include <cstring>
#include <getopt.h>

#include <algorithm>
#include <array>
#include <fstream>

#include "flavour.hpp"
#include "log.hpp"

namespace mc {
struct args {
  static void parse(int argc, char *argv[]) {
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
        log::info("Expecting a Max Clique of size", size);
        break;
      case 'u':
        upper_bound = static_cast<std::size_t>(std::atol(optarg));
        log::info("Expected Upper Bound for Max Clique of size", upper_bound);
        break;
      case 'l':
        lower_bound = static_cast<std::size_t>(std::atol(optarg));
        log::info("Expected Lower Bound for Max Clique of size", lower_bound);
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
      case 'h':
      case ':':
      case '?':
      default:
        std::cerr
            << "Find the Max Clique of a graph using branch-and-bound\n\n"
            << "  -r N run for N times\n"
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

    log::info("Running", args::exec_mode);

    if (stdin) {
      log::info("Reading from STDIN");
    } else {
      log::info("Reading from", filename);
    }
  }

  static inline std::istream &stream() { return stdin ? std::cin : file; }

  static inline std::ifstream file;
  static inline long num_turns = 5;
  static inline bool expect_size, undirected = true, stdin = true;
  static inline std::size_t size = -1u, upper_bound = -1u, lower_bound = 1;
  static inline std::string filename;
  static inline bool draw = false;

  static inline flavour exec_mode = flavour::heuristic;
};

struct input {
  inline input() {
    std::string line;
    std::getline(args::stream(), line);
    fetch_ftor fetcher(line);
    if (not fetcher(num_v) || not fetcher(num_e)) {
      throw std::runtime_error("Could not parse input header, use -?");
    }
    if (not args::expect_size && (args::expect_size = fetcher(args::size))) {
      log::info("Expecting a Max Clique of size", args::size);
    }

    log::info("Source Graph is", (args::undirected ? "Undirected" : "Directed"),
              "with", num_v, "vertices and", num_e, "edges");
  }

  ~input() { log::info("~input()"); }
  
  inline std::istream &operator*() { return args::stream(); }
  inline std::istream *operator->() { return &args::stream(); }

  inline bool is_maximal_size(std::size_t clique_size) const {
    return not args::expect_size || clique_size >= args::size;
  }

  struct fetch_ftor {
    explicit fetch_ftor(std::string &line) : iss(line) {}

    template <typename T> bool operator()(T &value) {
      if (iss.bad() || iss.eof()) {
        return false;
      } else {
        iss >> value;
        return true;
      }
    }

  private:
    std::istringstream iss;
  };

  std::size_t num_v, num_e;
};

struct feed {
  static inline const std::int64_t CHUNK_SIZE = 16384; // 16KB
  static inline const char deli = '\n', sep = ' ';

  using buffer = std::array<char, CHUNK_SIZE>;

  feed(feed &&feed) = delete;
  feed(const feed &feed) = delete;

  explicit inline feed(input &in) : in(in), tail_remaining(remaining.begin()) {
    std::memset(remaining.data(), 0x0, remaining.size());
  }

  ~feed() {
    if (tail_remaining != remaining.begin()) {
      log::error("INVALID file: no NL at the end of the file");
    }
  }

  inline operator bool() { return reading(); }
  inline std::size_t estimate_chunks() const {
    return 1 + in.num_e / CHUNK_SIZE;
  }
  inline std::size_t num_vertices() const { return in.num_v; }
  inline std::size_t num_edges() const { return in.num_e; }

  std::string read_chunk() {
    buffer buff;
    const buffer::difference_type size_remaining =
        std::distance(remaining.begin(), tail_remaining);
    std::move(remaining.begin(), tail_remaining, buff.data());

    auto read = read_next(CHUNK_SIZE - size_remaining);
    auto begin_read = read.first.begin();
    const auto size_read = read.second;
    std::move(begin_read, begin_read + size_read, buff.data() + size_remaining);

    const buffer::iterator tail_buff =
        buff.begin() + (size_remaining + size_read);
    buffer::iterator delimiter = tail_buff;
    while (delimiter != buff.begin() && *--delimiter != feed::deli)
      ;
    if (delimiter == buff.begin()) {
      log::error("FAILURE: chunk_size=", CHUNK_SIZE,
                 " exceeded! recompile with a bigger size");
    }

    const std::size_t buffer_size =
        static_cast<std::size_t>(std::distance(buff.begin(), delimiter++));
    std::move(delimiter, tail_buff, remaining.data());
    tail_remaining = remaining.begin() + std::distance(delimiter, tail_buff);

    return std::string{std::move(buff.data()), buffer_size};
  }

private:
  inline bool reading() {
    if (in->bad()) {
      log::error("UNEXPECTED READ FAILURE");
    }
    return not in->eof() && not in->fail();
  }

  inline std::pair<buffer, const std::streamsize>
  read_next(std::streamsize size_read) {
    buffer read;
    if (in->read(read.data(), size_read); in->bad()) {
      log::error("FAILURE: cannot read from stream");
    }
    return {read, in->gcount()};
  }

  input &in;
  buffer remaining{};
  buffer::iterator tail_remaining;
};
} // namespace mc

#endif // INPUT_HPP
