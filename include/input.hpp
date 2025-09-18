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

#include <atomic>
#include <cstdlib>
#include <cstring>
#include <getopt.h>

#include <array>
#include <filesystem>
#include <fstream>
#include <numeric>

#include "flavour.hpp"
#include "logger.hpp"

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

struct input_source {
  inline input_source() {
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

  ~input_source() { logger::debug("~input()"); }

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

#ifndef CHUNK_SIZE
#define CHUNK_SIZE (4 * 1024)
#endif

#define FEED_SCALE_JOBS 3

struct feed {
  static inline std::size_t feed_jobs_scale(std::size_t num_jobs) {
    return num_jobs * FEED_SCALE_JOBS;
  }

  static inline constexpr std::size_t BUFF_SIZE = CHUNK_SIZE * 1024;
  static_assert(BUFF_SIZE != 0);

  static inline constexpr char deli = '\n', sep = ' ';

  using buffer = std::pmr::vector<char>;

  feed(const feed &) = delete;
  feed &operator=(const feed &) = delete;

  explicit feed(input_source &in, std::size_t feed_size);

  ~feed();

  inline operator bool() { return reading(); }

  inline std::size_t num_vertices() const { return in.num_v; }
  inline std::size_t num_edges() const { return in.num_e; }

  using buffer_iterator = buffer::value_type *;

  struct chunk {
    inline void dispose() { owner->reclaim(idx); }

    inline std::size_t size() const {
      return static_cast<std::size_t>(end - begin);
    }

    const buffer_iterator begin;
    const buffer_iterator end;

  private:
    friend feed;

    inline chunk(buffer_iterator begin, buffer_iterator end, std::size_t idx,
                 feed *owner)
        : begin(begin), end(end), idx(idx), owner(owner) {}

    const std::size_t idx;
    feed *owner;
  };

  chunk read_chunk();

private:
  inline void reclaim(std::size_t idx) {
    std::scoped_lock lk(mtx);
    idxs.emplace_back(idx);
    logger::debug("reclaimed buffer", idx);
    recycle_cv.notify_one();
  }

  inline bool reading() const {
    if (in->bad()) {
      logger::error("UNEXPECTED READ FAILURE");
    }
    return !in->eof() && !in->fail();
  }

  input_source &in;

  buffer remaining{BUFF_SIZE, memory::pool()};
  buffer_iterator begin_rem = remaining.data();
  buffer_iterator tail_rem = remaining.data();

  std::vector<buffer> buffs;

  mutable std::mutex mtx;
  std::condition_variable recycle_cv;
  std::vector<std::size_t> idxs;
};
} // namespace mc

#endif // INPUT_HPP
