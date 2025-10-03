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

#include <array>
#include <filesystem>
#include <fstream>

#include "core.hpp"
#include "flavour.hpp"

namespace mc {
struct args {
  static void parse(int argc, char *argv[]);

  static inline std::istream &stream() {
    try {
      return stdin ? std::cin : file;
    } catch (std::runtime_error e) {
      std::cerr << "cannot open stream!\n";
      throw e;
    }
  }
  static inline std::size_t stream_size() {
    return stdin ? 1 : std::filesystem::file_size(filename);
  }

  static inline std::ifstream file;
  static inline long num_turns = 1;
  static inline bool expect_size, undirected = true, stdin = true;
  static inline std::size_t size = -1u, upper_bound = -1u, lower_bound = 1;
  static inline std::string filename;
  static inline bool draw = false;

  static std::uint16_t num_threads;
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
      logger::info("Expecting a Max Clique of size", args::size);
    }
    if (num_v == 0) {
      logger::error("Number of Vertices is Zero! abort.");
    }
    if (num_e == 0) {
      logger::error("Number of Edges is Zero! abort.");
    }
  }

  ~input() { logger::debug("~input()"); }

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
  static inline constexpr std::size_t CHUNK_SIZE = 2 * 1024 * 1024; // 1MB
  static_assert(CHUNK_SIZE != 0);

  static inline constexpr char deli = '\n', sep = ' ';

  using buffer = std::array<char, CHUNK_SIZE>;

  feed(feed &&feed) = delete;
  feed(const feed &feed) = delete;

  explicit inline feed(input &in)
      : in(in), begin_rem(remaining.data()), tail_rem(remaining.data()) {
    std::memset(remaining.data(), 0x00, remaining.size());
    const std::size_t _estimate_chunks = 1 + args::stream_size() / CHUNK_SIZE;
    const double avg_edge_line =
        static_cast<double>(args::stream_size()) / in.num_e;
    const std::size_t _edges_per_chunk = static_cast<std::size_t>(
        std::max(1.0, std::floor(CHUNK_SIZE / avg_edge_line)));
    logger::info("Stream Size", logger::size_unit(args::stream_size()),
                 "and Chunk Size", logger::size_unit(CHUNK_SIZE));
    logger::info("Estimating", logger::number_unit(_estimate_chunks),
                 "Chunks with", logger::number_unit(_edges_per_chunk),
                 "Edges per Chunk");
  }

  ~feed() {
    if (tail_rem != begin_rem) {
      // FIXME: use either exceptions or logger::error
      logger::error("INVALID file: no NL at the end of the file");
    }
    logger::debug("~feed()");
  }

  inline operator bool() { return reading(); }

  // inline std::size_t estimate_chunks() const { return _estimate_chunks; }
  // inline std::size_t estimate_num_edges() const { return _estimate_num_edges;
  // } inline std::size_t edges_per_chunk() const { return _edges_per_chunk; }
  inline std::size_t num_vertices() const { return in.num_v; }
  inline std::size_t num_edges() const { return in.num_e; }

  struct buffer_content {
    friend feed;

    buffer::const_iterator begin() const { return buff.begin(); }
    buffer::const_iterator end() const { return buff.begin() + size; }

  private:
    buffer buff;
    std::size_t size;
  };

  buffer_content read_chunk();

private:
  inline bool reading() {
    if (in->bad()) {
      logger::error("UNEXPECTED READ FAILURE");
    }
    return not in->eof() && not in->fail();
  }

  input &in;
  buffer remaining;
  char *begin_rem;
  char *tail_rem;

  // const std::size_t _estimate_chunks;
  // const std::size_t _estimate_num_edges;
  // const std::size_t _edges_per_chunk;
};
} // namespace mc
// namespace mc

#endif // INPUT_HPP
