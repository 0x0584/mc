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
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

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

#ifndef FEED_SCALE_JOBS
#define FEED_SCALE_JOBS 3
#endif

struct feed {
  static inline std::size_t feed_jobs_scale(std::size_t num_jobs) {
    return num_jobs * FEED_SCALE_JOBS;
  }

  static inline constexpr std::size_t BUFF_SIZE = CHUNK_SIZE * 1024;
  static_assert(BUFF_SIZE != 0);

  static inline constexpr char deli = '\n', sep = ' ';

  using buffer = char *;

  feed(const feed &) = delete;
  feed &operator=(const feed &) = delete;

  explicit feed(input_source &in, std::size_t num_jobs);

  ~feed();

  inline operator bool() {
    return !fetch_idxs.empty() || !done.load(std::memory_order_acquire);
  }

  inline std::size_t num_vertices() const { return in.num_v; }
  inline std::size_t num_edges() const { return in.num_e; }

  using buffer_iterator = buffer;

  struct chunk {
    inline void dispose() { owner->reclaim(idx); }

    inline std::size_t size() const {
      return static_cast<std::size_t>(end - begin);
    }

    buffer_iterator begin;
    buffer_iterator end;

    std::size_t idx;

    std::chrono::nanoseconds read_t;
    std::chrono::nanoseconds fetch_t;
    std::chrono::nanoseconds proc_t;
    std::chrono::nanoseconds dispos_t;

  private:
    friend feed;

    inline chunk(buffer_iterator begin, buffer_iterator end, std::size_t idx,
                 feed *owner)
        : begin(begin), end(end), idx(idx), owner(owner) {}

    feed *owner;
  };

  chunk fetch_chunk();

private:
  struct window {
    char *base = nullptr;
    size_t length = 0;
    off_t offset = 0;
    std::atomic_int refs = 0;
  };

  window map_window(int fd, off_t offset, size_t window_size,
                    size_t file_size) {
    if (static_cast<size_t>(offset) >= file_size) [[unlikely]] {
      return {nullptr, 0, offset};
    } else {
      file_size -= static_cast<size_t>(offset);
      const size_t len = std::min(window_size, file_size);
      void *addr = mmap(nullptr, len, PROT_READ, MAP_PRIVATE, fd, offset);
      if (addr == MAP_FAILED) [[unlikely]] {
        throw std::runtime_error("mmap failed");
      } else {
        madvise(addr, len, MADV_SEQUENTIAL);
        return {static_cast<char *>(addr), len, offset};
      }
    }
  }

  void unmap_window(window &w) {
    if (w.base) [[likely]] {
      munmap(w.base, w.length);
      w.base = nullptr;
      w.length = 0;
    }
  }

  void read_chunk(std::size_t idx);

  inline void reclaim(std::size_t idx) {
    thread_local auto last = std::chrono::high_resolution_clock::now();

    {
      std::scoped_lock lk(recycle_mtx);
      read_idxs.emplace_back(idx);
      auto now = std::chrono::high_resolution_clock::now();
      chunks[idx].dispos_t = now - last;
      last = now;
      logger::info("chunk", idx, "read", logger::duration(chunks[idx].read_t),
                   "fetch", logger::duration(chunks[idx].fetch_t), "disposed",
                   logger::duration(chunks[idx].dispos_t));
      recycle_cv.notify_one();
    }
  }

  inline bool reading() const {
    if (in->bad()) {
      logger::error("UNEXPECTED READ FAILURE");
    }
    return !in->eof() && !in->fail();
  }

  input_source &in;
  const std::size_t batch_size;
  const size_t win_read_size;
  const double prefetch_threshold = 0.75;

  buffer remaining;
  buffer_iterator begin_rem;
  buffer_iterator tail_rem;

  std::vector<buffer> buffs;
  std::vector<std::size_t> read_idxs;

  std::vector<chunk> chunks;
  std::vector<std::size_t> fetch_idxs;

  mutable std::mutex recycle_mtx;
  std::condition_variable recycle_cv;

  mutable std::mutex fetch_mtx;
  std::condition_variable fetch_cv;

  std::jthread reader;
  std::atomic_bool done = false;
};
} // namespace mc

#endif // INPUT_HPP
