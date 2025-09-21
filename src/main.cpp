// main.cpp
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

#include <chrono>
#include <climits>
#include <cstring>
#include <fcntl.h>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "graph.hpp"
#include "profiler.hpp"

// #include "mc.hpp"

class SystemInfo {
public:
  static inline size_t get_page_size() noexcept {
    static size_t page_size = calculate_page_size();
    return page_size;
  }

private:
  static inline size_t calculate_page_size() noexcept {
    if (long result = sysconf(_SC_PAGESIZE); result <= 0) [[unlikely]] {
      logger::warn("Failed to determine system page size via sysconf, fallback "
                   "4096 Bytes");
      return 4096u;
    } else {
      return static_cast<size_t>(result);
    }
  }
};

using namespace mc;

// === Config ===
// constexpr size_t CHUNK_SIZE   = 4 * 1024 * 1024;   // 4 MiB
// constexpr size_t WINDOW_SIZE = 256 * 1024 * 1024; // 256 MiB
// constexpr double PREFETCH_THRESHOLD = 0.75; // prefetch next window at 75%//
constexpr size_t INIT_WINDOW = 256 * 1024 * 1024; // 256 MiB
constexpr size_t MAX_WINDOW = 1024 * 1024 * 1024; // 1 GiB
constexpr double MIN_TH = 0.1;
constexpr double MAX_TH = 0.95;
constexpr double FAST_MS = 50.0;
constexpr double SLOW_MS = 200.0;

// === MMap Window ===
struct MMapWindow {
  char *base = nullptr;
  size_t length = 0;
  off_t offset = 0;
  std::future<void> prefetch;
};

MMapWindow map_window(int fd, off_t offset, size_t window_size,
                      size_t file_size) {
  if (offset >= (off_t)file_size)
    return {};
  size_t len = std::min(window_size, file_size - (size_t)offset);
  void *addr = mmap(nullptr, len, PROT_READ, MAP_PRIVATE, fd, offset);
  if (addr == MAP_FAILED)
    throw std::runtime_error("mmap failed");
  madvise(addr, len, MADV_WILLNEED);
  auto prefetch = std::async([addr, len]() {
    const char *ptr = static_cast<char *>(addr);
    for (size_t i = 0; i < len; i += SystemInfo::get_page_size()) {
      volatile char dummy = ptr[i];
      (void)dummy;
    }
  });
  return {static_cast<char *>(addr), len, offset, std::move(prefetch)};
}

void unmap_window(MMapWindow &w) {
  if (w.base) [[likely]] {
    if (w.prefetch.valid()) [[unlikely]] {
      w.prefetch.wait();
    }
    munmap(w.base, w.length);
    w = {};
  }
}

// === Reader ===
class Reader {
public:
  // === Chunk ===
  struct Chunk {
    inline Chunk() = default;
    inline Chunk(Reader *owner, size_t id) : owner(owner), id(id) {}

    inline void fill(const char *src, std::size_t sz) {
      if (size() < sz) [[unlikely]] {
        memory::deallocate<char>(first);
        first = memory::allocate<char>(sz);
      }
      std::memcpy(first, src, sz);
      last = first + sz;
    }

    inline void dispose() { owner->dispose(id); }
    inline size_t size() const { return static_cast<size_t>(last - first); }
    inline bool empty() const { return first == nullptr && last == nullptr; }
    inline const char *begin() const { return first; }
    inline const char *end() const { return last; }

  private:
    char *first = nullptr;
    char *last = nullptr;
    Reader *owner = nullptr;
    size_t id = static_cast<size_t>(-1);
  };

private:
  const std::size_t num_chunks;
  static const inline std::size_t num_hot_chunks = 8;

  std::mutex fetch_mtx;
  std::condition_variable fetch_cv;
  std::jthread fetch_thrd;

  std::pmr::vector<Chunk> chunks;
  std::pmr::vector<std::size_t> chunk_ids;

public:
  Reader(const std::string &filename, std::size_t num_jobs)
      : num_chunks(num_hot_chunks * num_jobs), chunks(memory::pool()),
        chunk_ids(memory::pool()) {
    fd_ = open(filename.c_str(), O_RDONLY);
    if (fd_ < 0) {
      throw std::runtime_error("open failed");
    }
    struct stat st;
    if (fstat(fd_, &st) < 0) {
      throw std::runtime_error("fstat failed");
    }
    file_size_ = static_cast<size_t>(st.st_size);

    chunks.reserve(num_chunks);
    chunk_ids.reserve(num_chunks);
    for (size_t i = 0; i < num_chunks; ++i) {
      chunks.emplace_back(this, i);
      chunk_ids.emplace_back(i);
    }
    window_size_ = INIT_WINDOW;
    current_ = map_window(fd_, 0, window_size_, file_size_);
    last_call_ = std::chrono::high_resolution_clock::now();
  }

  ~Reader() {
    logger::warn("total allocated chunks", chunks.size());
    unmap_window(current_);
    if (next_.valid()) {
      auto w = next_.get();
      unmap_window(w);
    }
    if (fd_ >= 0)
      close(fd_);
  }

  void dispose(size_t id) {
    {
      std::scoped_lock lk(fetch_mtx);
      chunk_ids.emplace_back(id);
    }
  }

  Chunk read_chunk() {
    auto now = std::chrono::high_resolution_clock::now();

    if (!current_.base) [[unlikely]] {
      return {};
    }

    if (processed_ >= current_.length) [[unlikely]] {
      unmap_window(current_);
      if (next_.valid()) {
        current_ = next_.get();
      }
      next_ = {};
      processed_ = 0;
      if (!current_.base) {
        return {};
      }
      logger::warn("used prefetched page");
    }

    auto start = std::chrono::high_resolution_clock::now();
    auto id_start = start, id_end = start;
    scope_dtor guard([&start, &id_start, &id_end] {
      auto end = std::chrono::high_resolution_clock::now();
      logger::info("chunk read in", logger::duration(start, end), "id took",
                   logger::duration(id_start, id_end));
    });

    size_t to_copy =
        std::min(static_cast<size_t>(CHUNK_SIZE),
                 static_cast<size_t>(current_.length - processed_));

    const char *src = current_.base + processed_;
    processed_ += to_copy;

    // double ms =
    //     std::chrono::duration<double, std::milli>(now - last_call_).count();

    id_start = std::chrono::high_resolution_clock::now();
    size_t chunk_id;
    {
      std::scoped_lock lk(fetch_mtx);
      if (chunk_ids.empty()) {
        chunk_id = chunks.size();
        chunk_ids.emplace_back(chunk_id);
        chunks.emplace_back(this, chunk_id);
      } else {
        chunk_id = chunk_ids.back();
        chunk_ids.pop_back();
      }
    }
    id_end = std::chrono::high_resolution_clock::now();

    // Copy out
    Chunk &chunk = chunks[chunk_id];
    chunk.fill(src, to_copy);

    double ms =
        std::chrono::duration<double, std::milli>(now - last_call_).count();
    last_call_ = now;

    // Adaptive threshold
    if (ms < FAST_MS) {
      threshold_ = std::min(MAX_TH, threshold_ * 1.1);
      fast_count_++;
    } else if (ms > SLOW_MS) {
      threshold_ = std::max(MIN_TH, threshold_ * 0.9);
      fast_count_ = 0;
    } else {
      fast_count_ = 0;
    }

    // Grow window if consistently fast
    if (fast_count_ >= 3 && window_size_ < MAX_WINDOW) {
      size_t new_size = std::min(MAX_WINDOW, (size_t)(window_size_ * 1.5));
      if (new_size > window_size_) {
        window_size_ = new_size;
        logger::warn("Increased window size to", window_size_);
      }
      fast_count_ = 0;
    }

    // Prefetch
    if (!next_.valid() && processed_ > (threshold_ * current_.length)) {
      off_t next_off = current_.offset + static_cast<off_t>(current_.length);
      if (next_off < (off_t)file_size_) {
        next_ = std::async(std::launch::async, [this, next_off] {
          return map_window(fd_, next_off, window_size_, file_size_);
        });
        logger::warn("Prefetched window at offset", next_off,
                     "threshold=", threshold_);
      }
    }

    return chunk;
  }

private:
  int fd_{-1};
  size_t file_size_{0};
  size_t window_size_{0};
  MMapWindow current_{};
  std::future<MMapWindow> next_;
  size_t processed_{0};

  double threshold_ = 0.5;
  int fast_count_ = 0;
  std::chrono::high_resolution_clock::time_point last_call_;
};

// === Simulated work ===
void process_chunk(const Reader::Chunk &c, size_t id) {
  auto start = std::chrono::high_resolution_clock::now();
  std::vector<char> v(c.begin(), c.end());
  auto end = std::chrono::high_resolution_clock::now();
  logger::info("Worker processed chunk", id, "with size",
               logger::size_unit(c.size()), "in", logger::duration(start, end));
}

int main(int argc, char *argv[]) {
  std::set_terminate([] {
    try {
      std::rethrow_exception(std::current_exception());
    } catch (const std::runtime_error &e) {
      logger::error(e.what());
    } catch (const std::exception &e) {
      logger::error(e.what());
    } catch (...) {
      logger::error("Uncaught unknown exception!");
    }
    std::abort(); // should never be reached!
  });

  args::parse(argc, argv);

  size_t total_chunks = 0;
  size_t total_bytes = 0;

  auto t0 = std::chrono::high_resolution_clock::now();

  // use our Reader abstraction
  Reader reader(args::filename, args::num_threads);

  // create pool with N threads
  thread::pool tasks(args::num_threads);

  profiler_start("mmap_read.prof");
  for (;;) {
    Reader::Chunk c = reader.read_chunk();
    if (c.empty())
      break;

    ++total_chunks;
    total_bytes += c.size();
    std::size_t chunk_id = total_chunks;

    // submit work to pool
    tasks.exec([c = std::move(c), chunk_id](auto) mutable {
      process_chunk(c, chunk_id); // or process_chunk(c) if you adapt signature
      c.dispose();
    });
  }
  tasks.join();
  profiler_stop();

  auto t1 = std::chrono::high_resolution_clock::now();

  logger::warn("Chunks", logger::size_unit(total_chunks), "total",
               logger::size_unit(total_bytes), "of avg.",
               logger::size_unit(CHUNK_SIZE), "Bytes", "Total time",
               logger::duration(t0, t1));

  return 42;

  /*
  std::ostringstream oss;
  oss << args::filename.c_str() << "_build.prof";
  std::string prof_file = oss.str();

  input_source in;
  graph_builder builder(in);
  profiler_start(prof_file.c_str());
  graph g = builder.build(args::undirected);
  profiler_stop();*/
  // g.print();

  // auto start = std::chrono::system_clock::now();

  // profiler_start("colour.prof");
  // auto cols = g.colour_sort(gen_keys(g.vertex_count()));
  // profiler_stop();

  // auto end = std::chrono::system_clock::now();
  // logger::warn("dsatur done in", logger::time_diff(start, end));

  // return 42;

  // multithreaded algo(g);
  // for (long turn = 1; turn <= args::num_turns; ++turn) {
  //   if (args::num_turns != 1) {
  //     logger::info("Turn", turn, "/", args::num_turns);
  //   }

  //   auto clique = algo.solve(args::exec_mode);
  //   std::pmr::vector<vertex> sorted_clique(clique.cbegin(), clique.cend(),
  //                                          memory::pool());
  //   std::sort(sorted_clique.begin(), sorted_clique.end());
  //   std::ostringstream oss;
  //   oss << "Max Clique has " << sorted_clique.size() << " vertices { ";
  //   for (vertex v : sorted_clique) {
  //     oss << v << " ";
  //   }
  //   oss << "}";
  //   logger::print(oss.str());

  //   if (args::draw) {
  //     algo.draw(clique);
  //   }
  // }

  return EXIT_SUCCESS;
}
