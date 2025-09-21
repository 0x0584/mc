#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstring>
#include <fcntl.h>
#include <functional>
#include <gperftools/profiler.h>
#include <iostream>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <utility>
#include <vector>

using namespace std::chrono;

// The purpose is to hide I/O latency, so we rely on the OS read-ahead.
// If the actual cost of I/O (the 'read') is hidden, the chunk processing
// becomes the bottleneck, which is the desired outcome for testing.

constexpr size_t CHUNK_SIZE = 1024 * 1024;
constexpr size_t WINDOW_SIZE = 128 * 1024 * 1024;
constexpr double PREFETCH_THRESHOLD = 0.10;

struct MMapWindow {
  char *base = nullptr;
  size_t length = 0;
  off_t offset = 0;
};

MMapWindow map_window(int fd, off_t offset, size_t window_size,
                      size_t file_size) {
  if (offset >= file_size)
    return {nullptr, 0, offset};

  size_t len = std::min(window_size, (size_t)(file_size - offset));

  void *addr = mmap(nullptr, len, PROT_READ, MAP_PRIVATE, fd, offset);
  if (addr == MAP_FAILED) {
    throw std::runtime_error("mmap failed");
  }

  int ret = madvise(addr, len, MADV_WILLNEED | MADV_SEQUENTIAL);
  if (ret != 0) {
    std::cerr << "Warning: madvise(MADV_WILLNEED) failed.\n";
  }

  return {static_cast<char *>(addr), len, offset};
}

void unmap_window(MMapWindow &w) {
  if (w.base) {
    madvise(w.base, w.length, MADV_DONTNEED);
    munmap(w.base, w.length);
    w.base = nullptr;
    w.length = 0;
  }
}

thread_local std::vector<char> buffer(CHUNK_SIZE);
void touch_chunk_memory(const char *first, const char *last) {
  std::memcpy(buffer.data(), first, last - first);
}

int main(int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " filename\n";
    return 1;
  }

  std::string filename = argv[1];
  int fd = open(filename.c_str(), O_RDONLY);
  if (fd < 0) {
    perror("open");
    return 1;
  }

  struct stat st;
  if (fstat(fd, &st) < 0) {
    perror("fstat");
    close(fd);
    return 1;
  }

  // ProfilerStart("mmap_read.prof");
  size_t file_size = st.st_size;

  size_t total_chunks = 0;
  size_t total_bytes = 0;

  auto t0 = std::chrono::high_resolution_clock::now();

  off_t offset = 0;

  MMapWindow current = map_window(fd, offset, WINDOW_SIZE, file_size);
  MMapWindow next = {nullptr, 0, 0};

  size_t processed_in_window = 0;
  while (current.base) {
    auto w_start = std::chrono::high_resolution_clock::now();

    while (processed_in_window < current.length) {
      auto c_start = std::chrono::high_resolution_clock::now();

      size_t to_copy =
          std::min(CHUNK_SIZE, current.length - processed_in_window);
      const char *first = current.base + processed_in_window;
      const char *last = first + to_copy;

      touch_chunk_memory(first, last);

      total_chunks++;
      total_bytes += to_copy;
      processed_in_window += to_copy;

      auto c_end = std::chrono::high_resolution_clock::now();
      auto c_dur =
          std::chrono::duration<double, std::milli>(c_end - c_start).count();
      std::cout << "    Chunk " << total_chunks << " took " << c_dur << " ms\n";

      if (!next.base &&
          processed_in_window > PREFETCH_THRESHOLD * current.length) {
        off_t next_off = current.offset + current.length;
        if (next_off < (off_t)file_size) {
          next = map_window(fd, next_off, WINDOW_SIZE, file_size);
          std::cout << "Prefetched window at offset " << next_off << "\n";
        }
      }
    }

    auto w_end = std::chrono::high_resolution_clock::now();
    auto w_dur =
        std::chrono::duration<double, std::milli>(w_end - w_start).count();
    std::cout << "Window at offset " << current.offset
              << " length=" << current.length << " processed in " << w_dur
              << " ms\n";

    unmap_window(current);
    current = next;
    next = {nullptr, 0, 0};
    processed_in_window = 0;
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  auto total_dur = std::chrono::duration<double, std::milli>(t1 - t0).count();

  std::cout << "\n=== Summary ===\n";
  std::cout << "File: " << filename << "\n";
  std::cout << "Size: " << file_size << " bytes\n";
  std::cout << "Chunks: " << total_chunks << " of ~" << CHUNK_SIZE
            << " bytes each\n";
  std::cout << "Processed: " << total_bytes << " bytes\n";
  std::cout << "Total time: " << total_dur << " ms\n";

  close(fd);

  // ProfilerStop();
  return 0;
}
