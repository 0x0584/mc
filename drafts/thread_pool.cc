
#define COL_GREEN "\x1b[32m"
#define COL_MAGENTA "\x1b[35m"
#define COL_CYAN "\x1b[36m"
#define COL_RESET "\x1b[0m"
#define COL_BOLD "\x1b[1m"

#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

#include <cassert>
#include <cstdint>

#include <functional>

#include <atomic>
#include <mutex>
#include <thread>

#include <queue>
#include <vector>

std::mutex print_mtx;

struct log {
  enum flags { none = 0x0, bold = 0b0001, ansi_colours = 0b0010 };

  static inline void setup_logger() {
    std::cout << std::fixed << std::setprecision(3) << std::left;
    std::cerr << std::fixed << std::setprecision(3) << std::left;
  }

  template <typename Chrono>
  static inline double duration(Chrono begin, Chrono end) {
    return std::chrono::duration<double>(end - begin).count();
  }

  template <typename Chrono>
  static inline std::string time_diff(Chrono begin, Chrono end,
                                      int flags = ansi_colours) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << std::left;
    if (flags & ansi_colours) {
      oss << COL_GREEN;
    }
    if (flags & bold) {
      oss << COL_BOLD;
    }
    oss << duration(begin, end) << "s";
    if (flags & ansi_colours) {
      oss << COL_RESET;
    }
    return oss.str();
  }

  static inline std::string progress(std::size_t index, std::size_t size) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << std::left
        << ((double(index + 1) * 100 / size)) << "%";
    return oss.str();
  }

  template <typename... Args> static inline void info(Args &&...args) {
    std::scoped_lock print_lock(mtx);
    std::cerr << "> ";
    ((std::cerr << std::forward<Args>(args) << " "), ...);
    std::cerr << COL_RESET << std::endl;
  }

  template <typename... Args> static inline void print(Args &&...args) {
    std::scoped_lock print_lock(mtx);
    ((std::cout << std::forward<Args>(args) << " "), ...);
    std::cout << COL_RESET << std::endl;
  }

  template <typename... Args> static inline void printv(Args &&...args) {
    std::scoped_lock print_lock(mtx);
    ((std::cout << std::forward<Args>(args)), ...); // verbose
    std::cout << COL_RESET << std::endl;
  }

  template <typename... Args> static inline void error(Args &&...args) {
    std::scoped_lock print_lock(mtx);
    ((std::cerr << std::forward<Args>(args) << " "), ...);
    std::cerr << COL_RESET << std::endl;
    std::exit(EXIT_FAILURE);
  }

  template <typename... Args>
  static inline void print_thread(std::uint32_t thread_id, Args &&...args) {
    std::scoped_lock print_lock(mtx);
    std::cout << std::string((thread_id + 1), ' ') << std::to_string(thread_id)
              << " ";
    ((std::cout << std::forward<Args>(args) << " "), ...);
    std::cout << COL_RESET << std::endl;
  }

private:
  static inline std::mutex mtx;
  // TODO: add a queue to handle incoming logs async
};

namespace thread {
const std::uint16_t num_threads = 8;

template <typename Task = std::function<void()>> struct pool {
  static inline const std::uint16_t max_num_threads = 255;

  pool(const pool &) = delete;
  pool(pool &&) = delete;

  explicit pool(std::uint16_t pool_size = 1) {
    assert(pool_size > 0);
    assert(pool_size <= std::min(max_num_threads, num_threads));
    _pool.resize(pool_size);
    for (std::uint16_t thread_id = 0u; thread_id < _pool.size(); ++thread_id) {
      _available.push(thread_id);
      { // FIXME: handle logs async
        std::unique_lock<std::mutex> lock(print_mtx);
        std::cerr << "thread " << thread_id << " is available\n";
      }
    }
    _joining = false;
  }

  ~pool() { join(); }

  pool &operator=(const pool &) = delete;
  pool &operator=(pool &&) = delete;

  void exec(Task &&task) {
    std::uint16_t thread_id;
    if (std::unique_lock<std::mutex> pool_lock(_pool_mtx);
        _joining || _available.empty()) {
      _pending_tasks.emplace(task);
      return;
    } else {
      thread_id = _available.front();
      _available.pop();
    }
    {
      std::unique_lock<std::mutex> print_lock(print_mtx);
      std::cerr << "thread " << thread_id << " is busy\n";
    }
    _pool[thread_id] =
        std::thread([this, thread_id, callback = std::move(task)]() {
          callback();
          _process_pending(thread_id);
          _available.emplace(thread_id);
          {
            std::unique_lock<std::mutex> print_lock(print_mtx);
            std::cerr << "thread " << thread_id << " is available\n";
          }
        });
  }

  void join(bool blocking = true) {
    _joining = true;
    {
      std::unique_lock<std::mutex> print_lock(print_mtx);
      std::cerr << "joining threads..\n";
    }
    std::thread barrier = std::thread([this]() {
      for (auto &task : _pool) {
        if (task.joinable()) {
          task.join();
        }
      }
      _joining = false;
      {
        std::unique_lock<std::mutex> print_lock(print_mtx);
        std::cerr << "all threads joined.\n";
      }
      assert(_pending_tasks.empty());
      assert(_available.size() == _pool.size());
    });

    if (blocking) {
      barrier.join();
    } else {
      barrier.detach();
      {
        std::unique_lock<std::mutex> print_lock(print_mtx);
        std::cerr << "non-blocking, detached.\n";
      }
    }
  }

private:
  void _process_pending(std::uint16_t thread_id) {
    while (true) {
      Task pending_task;
      if (std::unique_lock<std::mutex> pool_lock(_pool_mtx);
          _pending_tasks.empty()) {
        break;
      } else {
        pending_task = std::move(_pending_tasks.front());
        _pending_tasks.pop();
      }
      {
        std::unique_lock<std::mutex> lock(print_mtx);
        std::cerr << "thread " << thread_id
                  << " is processing a pending task\n";
      }
      pending_task();
    }
  }

  std::atomic_bool _joining;
  std::queue<Task> _pending_tasks;
  std::queue<std::uint16_t> _available;
  std::vector<std::thread> _pool;
  std::mutex _pool_mtx;
};
} // namespace thread

using namespace std::chrono;

struct Task {
  using duration_ms = std::chrono::duration<std::uint64_t, std::ratio<1, 1000>>;

  explicit Task(duration_ms d) : sleep_ms(d) { id = task_id++; }

  void operator()() {
    {
      std::unique_lock<std::mutex> print_lock(print_mtx);
      std::cerr << " #" << id << " began executing..\n";
    }

    std::this_thread::sleep_for(sleep_ms);

    {
      std::unique_lock<std::mutex> print_lock(print_mtx);
      std::cerr << " #" << id << " finished.\n";
    }
  }

private:
  static inline std::atomic_uint64_t task_id = 1;

  std::uint64_t id;
  duration_ms sleep_ms;
};

int main() {
  thread::pool pool(3);

  pool.exec(Task(19000ms));
  pool.exec(Task(1000ms));
  pool.exec(Task(100ms));
  pool.exec(Task(7000ms));
  pool.exec(Task(11000ms));

  pool.join();

  pool.exec(Task(900ms));
  pool.exec(Task(7000ms));
  pool.exec(Task(600ms));

  pool.join();

  pool.exec(Task(100ms));
  pool.exec(Task(500ms));
  pool.exec(Task(1000ms));
  pool.exec(Task(200ms));
  pool.exec(Task(400ms));
  pool.exec(Task(10000ms));
  pool.exec(Task(3000ms));

  pool.join();

  pool.exec(Task(31000ms));
  pool.exec(Task(5900ms));
  pool.exec(Task(800ms));
  pool.exec(Task(2000ms));
  pool.exec(Task(650ms));
  pool.exec(Task(100ms));
  pool.exec(Task(10ms));
  pool.exec(Task(20ms));
  pool.exec(Task(200ms));

  pool.join();

  pool.exec(Task(10000ms));
  pool.exec(Task(3000ms));
  pool.exec(Task(31000ms));
  pool.exec(Task(5900ms));
}
