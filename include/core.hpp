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

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <list>
#include <memory_resource>
#include <mutex>
#include <ostream>
#include <queue>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_set>

#include <cxxabi.h>
#include <execinfo.h>
#include <type_traits>

std::uint16_t __get_thread_id();

// since the thread can return on several conditions, it is
// practical to use a scope destructor to ensure that threads are
// marked as available not matter the branch
struct scope_dtor {
  scope_dtor(const scope_dtor &) = delete;
  scope_dtor(scope_dtor &&) = delete;

  inline explicit scope_dtor(std::function<void()> &&fn)
      : callback(std::move(fn)) {}
  inline ~scope_dtor() { callback(); }

  scope_dtor &operator=(const scope_dtor &) = delete;
  scope_dtor &operator=(scope_dtor &&) = delete;

private:
  std::function<void()> callback;
};

namespace memory {
std::pmr::synchronized_pool_resource *pool();

template <typename T, typename... Args>
inline std::shared_ptr<T> make_shared(Args &&...args) {
  std::pmr::polymorphic_allocator<T> alloc(pool());
  return std::allocate_shared<T>(alloc, std::forward<Args>(args)...);
}

template <typename T> struct deleter {
  inline deleter() noexcept : resource(nullptr) {}

  inline explicit deleter(std::pmr::memory_resource *res) noexcept
      : resource(res) {}

  inline void operator()(T *p) const {
    if (p && resource) {
      std::pmr::polymorphic_allocator<T> alloc(resource);
      alloc.deallocate(p, 1);
    }
  }

  inline const std::pmr::memory_resource *get_resource() const {
    return resource;
  }

private:
  std::pmr::memory_resource *resource;
};

template <typename T, typename... Args>
inline std::unique_ptr<T, deleter<T>> make_unique(Args &&...args) {
  void *memory = pool()->allocate(sizeof(T), alignof(T));
  T *object = new (memory) T(std::forward<Args>(args)...);
  return std::unique_ptr<T, deleter<T>>(object, PoolDeleter<T>(&pool));
}

struct gc {
  inline gc() : pool(options(), std::pmr::new_delete_resource()) {}

  inline std::pmr::polymorphic_allocator<std::byte> get_allocator() {
    return std::pmr::polymorphic_allocator<std::byte>(&pool);
  }

  template <typename T>
  inline std::pmr::polymorphic_allocator<T> get_allocator() {
    return std::pmr::polymorphic_allocator<T>(&pool);
  }

  inline std::pmr::synchronized_pool_resource *get_pool() { return &pool; }

private:
  static std::pmr::pool_options
  options(std::size_t max_blocks_per_chunk = 2,
          std::size_t largest_required_pool_block = 256) {
    std::pmr::pool_options opts;
    opts.max_blocks_per_chunk = max_blocks_per_chunk;
    opts.largest_required_pool_block = largest_required_pool_block;
    return opts;
  }

  std::pmr::synchronized_pool_resource pool;
};

} // namespace memory

template <typename T, typename = void> struct is_loggable : std::false_type {};

template <typename T>
struct is_loggable<T, std::void_t<decltype(std::declval<std::ostream &>()
                                           << std::declval<const T &>())>>
    : std::true_type {};

template <typename... Args> struct is_all_movable_helper : std::true_type {};

template <typename Head, typename... Tail>
struct is_all_movable_helper<Head, Tail...>
    : std::conditional_t<std::is_lvalue_reference_v<Head> ||
                             std::is_const_v<std::remove_reference_t<Head>>,
                         std::false_type, is_all_movable_helper<Tail...>> {};

template <typename... Args>
inline constexpr bool is_all_movable_v = is_all_movable_helper<Args...>::value;

struct logger {
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

  struct log_entry;

private:
  static inline std::mutex print_mtx;
  static inline std::mutex log_id_mtx;
  static inline std::size_t print_log_id = 1;

  static inline std::pmr::list<log_entry> log_queue;
  static inline std::mutex queue_mtx;
  static inline std::condition_variable logger_cv;
  static inline std::thread logger_thread;
  static inline bool stop_logger = false;

  static inline std::mutex flush_mtx;
  static inline std::condition_variable flush_cv;
  static inline std::atomic_bool flush_logs = false;

public:
  enum flags { none = 0x0, bold = 0b0001, ansi_colours = 0b0010 };

  static inline constexpr log_level current_level = log_level::LOG_LEVEL;

  static inline constexpr const char *get_level_str(log_level level) {
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

  template <typename Chrono>
  static inline auto duration(Chrono begin, Chrono end) {
    auto duration_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin)
            .count();
    std::ostringstream oss;
    if (duration_ns < 1000) {
      oss << duration_ns << " ns";
    } else if (duration_ns < 1000000) {
      double duration_us = static_cast<double>(duration_ns) / 1000.0;
      oss << duration_us << " µs";
    } else if (duration_ns < 1000000000) {
      double duration_ms = static_cast<double>(duration_ns) / 1000000.0;
      oss << duration_ms << " ms";
    } else {
      double duration_s = static_cast<double>(duration_ns) / 1000000000.0;
      oss << duration_s << " s";
    }
    return oss.str();
  }

  template <typename Chrono>
  static inline std::string time_diff(Chrono begin, Chrono end,
                                      int flags = ansi_colours) {
    (void)flags;
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << std::left;
    oss << duration(begin, end);
    return oss.str();
  }

  static inline std::string progress(std::size_t index, std::size_t size) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << std::left
        << ((double(index + 1) * 100 / size)) << "%";
    return oss.str();
  }

  static void flush() {
    {
      std::scoped_lock queue_lock(queue_mtx);
      flush_logs = true;
      logger_cv.notify_one();
    }
    {
      std::unique_lock flush_lock(flush_mtx);
      flush_cv.wait(flush_lock, [] { return !flush_logs; });
    }
  }

  static std::string stacktrace() {
    const int max_frames = 128;
    std::vector<void *> callstack(max_frames);
    int frames = backtrace(callstack.data(), max_frames);
    char **symbols = backtrace_symbols(callstack.data(), frames);

    std::ostringstream oss;
    oss << "\n";
    if (symbols) {
      for (int i = 0; i < frames; ++i) {
        char *demangled_name = nullptr;
        size_t len = 0;
        int status;
        char *mangled_name = symbols[i];

        demangled_name =
            abi::__cxa_demangle(mangled_name, demangled_name, &len, &status);

        if (status == 0 && demangled_name) {
          oss << demangled_name << '\n';
          free(demangled_name);
        } else {
          oss << symbols[i] << '\n';
        }
      }
      free(symbols);
    }
    return oss.str();
  }

  template <typename... Args>
  static inline std::string process_print_message(Args &&...args) {
    std::ostringstream oss;
    ((oss << std::forward<Args>(args) << " "), ...);
    oss << '\n';
    return oss.str();
  }

  template <typename... Args>
  static inline std::string process_printv_message(Args &&...args) {
    std::ostringstream oss;
    ((oss << std::forward<Args>(args)), ...);
    oss << '\n';
    return oss.str();
  }

  template <typename... Args>
  static inline std::string process_message(auto Level, auto now,
                                            auto colour_code, auto thread_id,
                                            Args &&...args) {
    std::ostringstream oss;
    oss << colour_code;
    oss << now << " [" << std::setfill('0') << std::setw(3) << thread_id << "] "
        << get_level_str(Level) << " ";
    ((oss << std::forward<Args>(args) << " "), ...);
    oss << COL_RESET << '\n';
    return oss.str();
  }

  template <log_level Level, typename... Args>
  static inline void _log_impl(const char *colour_code, Args &&...log_args) {
    static_assert((is_loggable<std::decay_t<Args>>::value && ...),
                  "operator<< overload missing.");
    if constexpr (Level <= current_level) {
      auto now = std::chrono::system_clock::now();

      std::function<std::string()> message;
      if constexpr (is_all_movable_v<Args...>) {
        message = [now, colour_code, thread_id = __get_thread_id(),
                   ... log_args_captured = std::forward<Args>(log_args)] {
          return process_message(Level, now, colour_code, thread_id,
                                 log_args_captured...);
        };
      } else {
        std::string processed_msg =
            process_message(Level, now, colour_code, __get_thread_id(),
                            std::forward<Args>(log_args)...);
        message = [msg = std::move(processed_msg)] { return msg; };
      }

      std::scoped_lock lock(queue_mtx);
      log_queue.emplace_back(std::move(message));
    }
  }

  [[nodiscard]] static scope_dtor setup_logger();

  template <typename... Args> static void debug(Args &&...args) {
    _log_impl<log_level::debug>(COL_CYAN, std::forward<Args>(args)...);
  }

  template <typename... Args> static void info(Args &&...args) {
    _log_impl<log_level::info>(COL_RESET, std::forward<Args>(args)...);
  }

  template <typename... Args> static void warn(Args &&...args) {
    _log_impl<log_level::warn>(COL_YELLOW, std::forward<Args>(args)...);
  }

  template <typename... Args> static void error(Args &&...args) {
    _log_impl<log_level::error>(COL_RED, std::forward<Args>(args)...,
                                stacktrace(), '\n');
    std::exit(EXIT_FAILURE);
  }

  template <typename... Args> static void print(Args &&...log_args) {
    static_assert((is_loggable<std::decay_t<Args>>::value && ...),
                  "operator<< overload missing.");

    std::function<std::string()> message;
    if constexpr (is_all_movable_v<Args...>) {
      message = [... log_args_captured = std::forward<Args>(log_args)] {
        return process_print_message(log_args_captured...);
      };
    } else {
      std::string processed_msg =
          process_print_message(std::forward<Args>(log_args)...);
      message = [msg = std::move(processed_msg)] { return msg; };
    }

    std::scoped_lock lock(queue_mtx);
    log_queue.emplace_back(std::move(message));
  }

  template <typename... Args> static void printv(Args &&...log_args) {
    static_assert((is_loggable<std::decay_t<Args>>::value && ...),
                  "operator<< overload missing.");

    std::function<std::string()> message;
    if constexpr (is_all_movable_v<Args...>) {
      message = [... log_args_captured = std::forward<Args>(log_args)] {
        return process_print_message(log_args_captured...);
      };
    } else {
      std::string processed_msg =
          process_print_message(std::forward<Args>(log_args)...);
      message = [msg = std::move(processed_msg)] { return msg; };
    }

    std::scoped_lock lock(queue_mtx);
    log_queue.emplace_back(std::move(message));
  }

  static inline const auto logger_destoy = logger::setup_logger();
}; // namespace logger

struct logger::log_entry {
  friend std::less<log_entry>;
  friend std::greater<log_entry>;
  friend std::equal_to<log_entry>;

  log_entry() = default;

  log_entry(const log_entry &) = delete;
  log_entry(log_entry &&) = default;

  log_entry &operator=(const log_entry &) = delete;
  log_entry &operator=(log_entry &&) = default;

  log_entry(std::function<std::string()> &&msg) : msg(std::move(msg)) {
    log_id = ++global_log_id;
  }

  std::size_t id() const { return log_id; }
  std::string message() const { return msg(); }

private:
  static inline std::atomic_size_t global_log_id = 0;

  std::size_t log_id; // XXX: add id_to_print
  std::function<std::string()> msg;
};

namespace std {
template <> struct less<logger::log_entry> {
  bool operator()(const logger::log_entry &lhs, const logger::log_entry &rhs) {
    return lhs.log_id < rhs.log_id;
  }
};

template <> struct greater<logger::log_entry> {
  bool operator()(const logger::log_entry &lhs, const logger::log_entry &rhs) {
    return lhs.log_id > rhs.log_id;
  }
};

template <> struct equal_to<logger::log_entry> {
  bool operator()(const logger::log_entry &lhs, const logger::log_entry &rhs) {
    return lhs.log_id == rhs.log_id;
  }
};
} // namespace std

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

namespace thread {
extern const std::uint16_t threads_per_core;
extern const std::uint16_t num_available_threads;

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

struct scope_timer {
  explicit inline scope_timer(const char *msg)
      : callback([msg = std::move(msg),
                  start = std::chrono::high_resolution_clock::now()] {
          auto end = std::chrono::high_resolution_clock::now();
          logger::warn(std::move(msg), logger::time_diff(start, end));
        }) {}

private:
  scope_dtor callback;
};

// #ifndef NDEBUG
#define make_scope_timer(x) scope_timer x(#x)
// #else
// #define make_scope_timer(x) \
//   do { \ } while (0)
// #endif

// FIXME: turn logger into a class

#endif // CORE_HPP
