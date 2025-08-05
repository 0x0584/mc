// core.hpp
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

#ifndef CORE_HPP
#define CORE_HPP

#define COL_RESET "\x1b[0m"
#define COL_BOLD "\x1b[1m"
#define COL_FAINT "\x1b[2m"
#define COL_ITALIC "\x1b[3m"
#define COL_UNDERLINE "\x1b[4m"

#define COL_BLACK "\x1b[30m"
#define COL_RED "\x1b[31m"
#define COL_GREEN "\x1b[32m"
#define COL_YELLOW "\x1b[33m"
#define COL_BLUE "\x1b[34m"
#define COL_MAGENTA "\x1b[35m"
#define COL_CYAN "\x1b[36m"
#define COL_WHITE "\x1b[37m"
#define COL_DEFAULT_FG "\x1b[39m"

#define BG_BLACK "\x1b[40m"
#define BG_RED "\x1b[41m"
#define BG_GREEN "\x1b[42m"
#define BG_YELLOW "\x1b[43m"
#define BG_BLUE "\x1b[44m"
#define BG_MAGENTA "\x1b[45m"
#define BG_CYAN "\x1b[46m"
#define BG_WHITE "\x1b[47m"
#define BG_DEFAULT_BG "\x1b[49m"

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <ostream>
#include <queue>
#include <sstream>
#include <string>
#include <thread>

namespace {
std::atomic_uint16_t __thread_id_count = 0;
thread_local bool __thread_id_assigned = false;
thread_local std::uint16_t __thread_id = 0;
} // namespace

inline std::uint16_t __get_thread_id() {
  if (not __thread_id_assigned) {
    __thread_id = __thread_id_count.fetch_add(1, std::memory_order_relaxed);
    __thread_id_assigned = true;
  }
  return __thread_id;
}

// since the thread can return on several conditions, it is
// practical to use a scope destructor to ensure that threads are
// marked as available not matter the branch
template <typename Callable> struct scope_dtor {
  scope_dtor(const scope_dtor &) = delete;
  scope_dtor(scope_dtor &&) = delete;

  inline explicit scope_dtor(Callable &&fn)
      : callback(std::forward<Callable>(fn)) {}
  inline ~scope_dtor() { callback(); }

  scope_dtor &operator=(const scope_dtor &) = delete;
  scope_dtor &operator=(scope_dtor &&) = delete;

private:
  Callable callback;
};

namespace logger {
enum class log_level : unsigned {
  off = 0,
  error,
  warn,
  info,
  debug

#ifndef LOG_LEVEL
#define LOG_LEVEL info
#endif
};

template <typename... Args> inline void debug(Args &&...args);
template <typename... Args> inline void info(Args &&...args);
template <typename... Args> inline void warn(Args &&...args);
template <typename... Args> inline void error(Args &&...args);
template <typename... Args> inline void print(Args &&...args);
template <typename... Args> inline void printv(Args &&...args);

#ifdef assert
#undef assert
#endif
#ifdef NDEBUG
#define assert(cond, ...) ((void)0)
#else
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define assert(cond, ...)                                                      \
  ((void)(cond ? ((void)0)                                                     \
               : logger::error(#cond, "failed, hint:", ##__VA_ARGS__, "in",    \
                               __func__, "at",                                 \
                               __FILE__ ":" TOSTRING(__LINE__))))
#endif

template <log_level level = log_level::info, typename... Args>
inline void print_thread(std::uint16_t thread_id, const char *colour_code,
                         Args &&...args);
} // namespace logger

namespace thread {
#ifndef THREADS_PER_CORE
#define THREADS_PER_CORE 2
#endif
constexpr std::uint16_t threads_per_core = THREADS_PER_CORE;
const std::uint16_t num_available_threads =
    std::thread::hardware_concurrency() * threads_per_core;

struct pool {
  using Task = std::function<void(std::uint16_t)>;

  pool(const pool &) = delete;
  pool(pool &&) = delete;

  explicit pool(std::uint16_t pool_size) {
    assert(pool_size > 0, "pool size cannot be 0");
    assert(pool_size <= thread::num_available_threads, "pool size", pool_size,
           "is greater than available threads", num_available_threads);
    _pool.resize(pool_size);
    for (std::uint16_t task_id = 0; task_id < _pool.size(); ++task_id) {
      _available.push(task_id);
    }
  }

  static pool sole_thread() { return pool(1); }

  ~pool() { logger::debug("~pool()"); }

  pool &operator=(const pool &) = delete;
  pool &operator=(pool &&) = delete;

  void exec(Task &&task) {
    std::uint16_t task_id;
    if (std::unique_lock<std::mutex> pool_lock(_pool_mtx); _available.empty()) {
      _pending_tasks.emplace(std::move(task));
      pool_lock.unlock();
      if (not pool_is_full) { // FIXME: use timed condition value to wait
        logger::debug("pool is full!");
        pool_is_full = true;
      }
      return;
    } else {
      task_id = _available.front();
      _available.pop();
      pool_is_full = false;
    }
    _pool[task_id] = std::jthread(
        [this, task_id, callback = std::forward<Task>(task)] mutable {
          callback(task_id);
          _process_pending(task_id);
        });
  }

  void join() {
    logger::debug("joining threads..");
    for (auto &task : _pool) {
      if (task.joinable()) {
        task.join();
      }
    }
    logger::debug("all threads joined.");
  }

  [[nodiscard]] bool pool_full() const { return pool_is_full; }

  void discard_pending() {
    std::unique_lock pool_lock(_pool_mtx);
    _pending_tasks = std::queue<Task>();
  }

private:
  inline void _process_pending(std::uint16_t task_id) {
    while (true) {
      std::unique_lock<std::mutex> pool_lock(_pool_mtx);
      if (_pending_tasks.empty()) { // FIXME: check for data race
        break;
      }
      Task pending_task = std::move(_pending_tasks.front());
      _pending_tasks.pop();
      pool_lock.unlock();
      pending_task(task_id);
    }
    _available.emplace(task_id);
  }

  std::queue<Task> _pending_tasks;
  std::queue<std::uint16_t> _available;
  std::vector<std::jthread> _pool;
  std::mutex _pool_mtx;
  bool pool_is_full = false;
};
} // namespace thread

namespace logger {
struct log_entry {
  using time_point = std::chrono::time_point<std::chrono::system_clock>;

  friend std::less<log_entry>;
  friend std::greater<log_entry>;
  friend std::equal_to<log_entry>;

  log_entry() = default;
  log_entry(time_point stamp, std::future<std::string> msg)
      : stamp(stamp), msg(std::move(msg)) {}

  std::string message() const { return msg.get(); }
  const time_point &timestamp() const { return stamp; }

private:
  time_point stamp;
  mutable std::future<std::string> msg;
};
} // namespace logger

namespace std {
template <> struct less<logger::log_entry> {
  bool operator()(const logger::log_entry &lhs, const logger::log_entry &rhs) {
    return lhs.stamp < rhs.stamp;
  }
};

template <> struct greater<logger::log_entry> {
  bool operator()(const logger::log_entry &lhs, const logger::log_entry &rhs) {
    return lhs.stamp > rhs.stamp;
  }
};

template <> struct equal_to<logger::log_entry> {
  bool operator()(const logger::log_entry &lhs, const logger::log_entry &rhs) {
    return lhs.stamp == rhs.stamp;
  }
};
} // namespace std

namespace logger {
static inline std::mutex print_mtx;
using log_queue_t = std::priority_queue<log_entry, std::pmr::vector<log_entry>,
                                        std::greater<log_entry>>;
static inline log_queue_t log_queue;
static inline std::mutex queue_mtx;
static inline std::condition_variable logger_cv;
static inline std::thread logger_thread;
static inline bool stop_logger = false;

enum flags { none = 0x0, bold = 0b0001, ansi_colours = 0b0010 };

static constexpr log_level current_level = log_level::LOG_LEVEL;

inline constexpr const char *get_level_str(log_level level) {
  switch (level) {
  case log_level::debug:
    return "DEBUG";
  case log_level::info:
    return "INFO ";
  case log_level::warn:
    return "WARN ";
  case log_level::error:
    return "ERROR";
  default:
    throw std::logic_error("unknown log_level!");
  }
}

[[nodiscard]] static inline auto setup_logger() {
  std::cout << std::fixed << std::setprecision(3) << std::left;
  std::cerr << std::fixed << std::setprecision(3) << std::left;
  logger_thread = std::thread([] {
    while (true) {
      std::unique_lock queue_lock(queue_mtx);
      logger_cv.wait(queue_lock, [] {
        std::unique_lock print_lock(print_mtx);
        return stop_logger || not log_queue.empty();
      });

      if (stop_logger && log_queue.empty()) {
        break;
      }

      std::string message = log_queue.top().message();
      log_queue.pop();
      queue_lock.unlock();

      std::unique_lock print_lock(print_mtx);
      std::cerr << message;
    }
  });
  return scope_dtor([]() {
    {
      std::unique_lock<std::mutex> lock(queue_mtx);
      stop_logger = true;
    }
    logger_cv.notify_one();
    if (logger_thread.joinable()) {
      logger_thread.join();
    }
  });
}

const auto logger_destoy = logger::setup_logger();

template <typename Chrono> inline double duration(Chrono begin, Chrono end) {
  return std::chrono::duration<double>(end - begin).count();
}

template <typename Chrono>
inline std::string time_diff(Chrono begin, Chrono end,
                             int flags = ansi_colours) {
  (void)flags;
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(3) << std::left;
  // if (flags & ansi_colours) {
  //   oss << COL_GREEN;
  // }
  // if (flags & bold) {
  //   oss << COL_BOLD;
  // }
  oss << duration(begin, end) << "s";
  // if (flags & ansi_colours || flags & bold) {
  //   oss << COL_RESET;
  // }
  return oss.str();
}

inline std::string progress(std::size_t index, std::size_t size) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2) << std::left
      << ((double(index + 1) * 100 / size)) << "%";
  return oss.str();
}

template <typename... Args>
struct is_copy_constructible_helper : std::true_type {};

template <typename Head, typename... Tail>
struct is_copy_constructible_helper<Head, Tail...>
    : std::conditional_t<std::is_copy_constructible_v<std::decay_t<Head>>,
                         is_copy_constructible_helper<Tail...>,
                         std::false_type> {};

template <typename... Args>
inline constexpr bool is_copy_constructible_v =
    is_copy_constructible_helper<Args...>::value;

template <log_level Level, typename... Args>
static inline void _log_impl(const char *colour_code, Args &&...log_args) {
  if constexpr (Level <= current_level) {
    auto now = std::chrono::system_clock::now();
    auto log_message = [now, colour_code,
                        thread_id = __get_thread_id()](auto &&...args) {
      auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                    now.time_since_epoch()) %
                1000;
      std::time_t tt = std::chrono::system_clock::to_time_t(now);
      std::tm tm_buf;
      localtime_r(&tt, &tm_buf);

      std::ostringstream oss;
      oss << colour_code;
      oss << std::put_time(&tm_buf, "%Y-%m-%d %H:%M:%S") << '.'
          << std::setfill('0') << std::setw(3) << ms.count() << " ["
          << std::setfill('0') << std::setw(3) << thread_id << "] "
          << get_level_str(Level) << " ";

      ((oss << std::forward<decltype(args)>(args) << " "), ...);

      oss << COL_RESET << '\n';
      return oss.str();
    };

    std::unique_lock lock(queue_mtx);
    if constexpr (is_copy_constructible_v<Args...>) {
      log_queue.emplace(now, std::async(std::launch::deferred, log_message,
                                        std::forward<Args>(log_args)...));
    } else {
      log_queue.emplace(
          now, std::async(
                   std::launch::async,
                   [](std::string message) { return std::move(message); },
                   log_message(std::forward<Args>(log_args)...)));
    }
    logger_cv.notify_one();
  }
}

template <typename... Args> inline void debug(Args &&...args) {
  _log_impl<log_level::debug>(COL_CYAN, std::forward<Args>(args)...);
}

template <typename... Args> inline void info(Args &&...args) {
  _log_impl<log_level::info>(COL_RESET, std::forward<Args>(args)...);
}

template <typename... Args> inline void warn(Args &&...args) {
  _log_impl<log_level::warn>(COL_YELLOW, std::forward<Args>(args)...);
}

template <typename... Args> inline void error(Args &&...args) {
  _log_impl<log_level::error>(COL_RED, std::forward<Args>(args)...);
  std::exit(EXIT_FAILURE);
}

template <typename... Args> inline void print(Args &&...args) {
  std::scoped_lock print_lock(print_mtx);
  ((std::cout << std::forward<Args>(args) << " "), ...);
  std::cout << COL_RESET << std::endl;
}

template <typename... Args> inline void printv(Args &&...args) {
  std::scoped_lock print_lock(print_mtx);
  ((std::cout << std::forward<Args>(args)), ...); // verbose
  std::cout << COL_RESET << std::endl;
}

template <log_level level, typename... Args>
inline void print_thread(std::uint32_t thread_id, const char *colour_code,
                         Args &&...args) {
  if constexpr (level <= current_level) {
    std::scoped_lock print_lock(print_mtx);
    std::cout << std::string((thread_id + 1), ' ') << colour_code
              << std::setw(2) << thread_id << COL_RESET << " ";
    ((std::cout << std::forward<Args>(args) << " "), ...);
    std::cout << COL_RESET << std::endl;
  }
}

} // namespace logger

#endif // CORE_HPP
