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

// TODO: create utils namespace

// since the thread can return on several conditions, it is
// practical to use a scope destructor to ensure that threads are
// marked as available not matter the branch
// FIXME: rename to scope_guard
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
// TODO: add a variant for std::pmr::unsynchronized_pool_resource
// std::pmr::unsynchronized_pool_resource *pool_unsafe();
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

// FIXME: rename this to pool
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

namespace std {
template <typename T, typename U>
inline ostream &operator<<(ostream &oss, const pair<T, U> &p) {
  return oss << "{" << p.first << ", " << p.second << "}";
}
} // namespace std

// TODO: unify the coding style
struct logger {
  struct colours {
    struct attr {
      static constexpr const char *reset = "\x1b[0m";
      static constexpr const char *bold = "\x1b[1m";
      static constexpr const char *faint = "\x1b[2m";
      static constexpr const char *italic = "\x1b[3m";
      static constexpr const char *underline = "\x1b[4m";
    };

    struct fg {
      static constexpr const char *black = "\x1b[30m";
      static constexpr const char *red = "\x1b[31m";
      static constexpr const char *green = "\x1b[32m";
      static constexpr const char *yellow = "\x1b[33m";
      static constexpr const char *blue = "\x1b[34m";
      static constexpr const char *magenta = "\x1b[35m";
      static constexpr const char *cyan = "\x1b[36m";
      static constexpr const char *white = "\x1b[37m";
      static constexpr const char *standard = "\x1b[39m";
    };

    struct bg {
      static constexpr const char *black = "\x1b[40m";
      static constexpr const char *red = "\x1b[41m";
      static constexpr const char *green = "\x1b[42m";
      static constexpr const char *yellow = "\x1b[43m";
      static constexpr const char *blue = "\x1b[44m";
      static constexpr const char *magenta = "\x1b[45m";
      static constexpr const char *cyan = "\x1b[46m";
      static constexpr const char *white = "\x1b[47m";
      static constexpr const char *standard = "\x1b[49m";
    };

    template <typename... Args>
    static inline std::string apply(Args &&...args) {
      std::ostringstream oss;
      ((oss << std::forward<Args>(args)), ...);
      oss << attr::reset;
      return oss.str();
    }
  };

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
  static inline std::mutex stderr_mtx;
  static inline std::mutex stdout_mtx;
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
      return "?";
    }
  }

  template <typename T> static inline std::string trim_trailing_zeros(T val) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << std::left << val;
    std::string fmt_val = oss.str();
    if constexpr (std::is_floating_point<T>::value) {
      while (fmt_val.back() == '0' && fmt_val.size() > 1) {
        fmt_val.pop_back();
      }
      if (fmt_val.back() == '.') {
        fmt_val.pop_back();
      }
    }
    return fmt_val;
  }

  template <typename Chrono>
  static inline auto duration(Chrono begin, Chrono end) {
    static constexpr std::array<std::pair<long long, const char *>, 7> units = {
        {{86'400'000'000'000LL, "day"},
         {3'600'000'000'000LL, " hours"},
         {60'000'000'000LL, " minutes"},
         {1'000'000'000LL, " s"},
         {1'000'000LL, " ms"},
         {1000LL, " µs"},
         {1LL, " ns"}}};
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << std::left;
    auto stamp =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin)
            .count();
#pragma unroll 8
    for (unsigned i = 0; i < units.size(); ++i) {
      if (stamp >= units[i].first) {
        double value = static_cast<double>(stamp) / units[i].first;
        oss << value << " " << units[i].second << (value < 2.0 ? "" : "s");
        break;
      }
    }
    return oss.str();
  }

  template <typename Units, typename Value>
  static inline std::string scale(Value val, double scale_factor,
                                  const Units &units) {
    double d_val = static_cast<double>(val);
    unsigned unit_index = 0;
#pragma unroll 8
    while (val >= scale_factor && unit_index < units.size() - 1) {
      val /= scale_factor;
      unit_index++;
    }
    std::ostringstream oss;
    oss << trim_trailing_zeros(d_val) << " " << units[unit_index];
    return oss.str();
  }

  static inline std::string number_unit(std::size_t val) {
    if (val < 1000)
      return trim_trailing_zeros(val);
    static constexpr std::array<char, 5> units = {'?', 'K', 'M', 'B', 'T'};
    static constexpr double scale_factor = 1000.0;
    return scale(val, scale_factor, units);
  }

  static inline std::string size_unit(std::size_t val) {
    static constexpr std::array<std::string, 5> units = {"B", "KB", "MB", "GB",
                                                         "TB"};
    static constexpr double scale_factor = 1024.0;
    return scale(val, scale_factor, units);
  }

  static inline std::string throughput(double throu) {
    static constexpr std::array<std::string, 5> units = {"/s", "K/s", "M/s",
                                                         "B/s", "T/s"};
    static constexpr double scale_factor = 1000.0;
    return scale(throu, scale_factor, units);
  }

  template <typename SizeType, typename Chrono>
  static inline std::string throughput(Chrono begin, Chrono end,
                                       SizeType count) {
    const double thrpt =
        count / std::chrono::duration<double>(end - begin).count();
    return throughput(thrpt);
  }

  // FIXME: remove this, deprecated! duration

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
        << trim_trailing_zeros((double(index + 1) * 100 / size)) << "%";
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
  static inline std::string process_message(auto Level, auto now,
                                            auto colour_code, auto thread_id,
                                            Args &&...args) {
    std::ostringstream oss;
    oss << colour_code;
    oss << now << " [" << std::setfill('0') << std::setw(3) << thread_id << "] "
        << get_level_str(Level) << " ";
    ((oss << std::forward<Args>(args) << " "), ...);
    oss << colours::attr::reset << '\n';
    return oss.str();
  }

  static inline void puts(const std::string &s, auto strm) {
    puts(s.c_str(), strm);
  }

  static inline void puts(const char *s, auto strm) { std::fputs(s, strm); }

  // non-blocking log for trivially movable types, otherwise, they are processed
  // in-place to avoid the overhead of potential expensive copy.  logs are sent
  // to `stderr'.
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
    _log_impl<log_level::debug>(colours::fg::cyan, std::forward<Args>(args)...);
  }

  template <typename... Args> static void info(Args &&...args) {
    _log_impl<log_level::info>(colours::attr::reset,
                               std::forward<Args>(args)...);
  }

  template <typename... Args> static void warn(Args &&...args) {
    _log_impl<log_level::warn>(colours::fg::yellow,
                               std::forward<Args>(args)...);
  }

  template <typename... Args> static void error(Args &&...args) {
    _log_impl<log_level::error>(colours::fg::red, std::forward<Args>(args)...,
                                stacktrace(), '\n');
    std::exit(EXIT_FAILURE);
  }

  // prints are processed in-place and displayed immediately to `stdout'
  // separated by spaces
  template <typename... Args> static void print(Args &&...log_args) {
    static_assert((is_loggable<std::decay_t<Args>>::value && ...),
                  "operator<< overload missing.");
    std::scoped_lock lk(stdout_mtx);
    puts(process_print_message(std::forward<Args>(log_args)...), stdout);
  }

  // verbose print is similar to `print()', but does not separate the parameters
  // by spaces.
  template <typename... Args> static void printv(Args &&...log_args) {
    static_assert((is_loggable<std::decay_t<Args>>::value && ...),
                  "operator<< overload missing.");
    std::scoped_lock lk(stdout_mtx);
    puts(colours::apply(std::forward<Args>(log_args)..., '\n'), stdout);
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

#ifdef NDEBUG
#define Assert(cond, ...) ((void)0)
#else
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define Assert(cond, ...)                                                      \
  ((void)(cond ? ((void)0)                                                     \
               : logger::error(#cond, "failed, hint:", ##__VA_ARGS__, "in",    \
                               __func__, "at",                                 \
                               __FILE__ ":" TOSTRING(__LINE__))))
#endif

namespace thread {
extern const std::uint16_t threads_per_core;
extern const std::uint16_t num_available_threads;

// fixed-size thread pool. it accepts tasks via `exec()' which are executed
// immediately after they are submitted. if a task is submitted while all
// threads are busy, it is enqueued and immediately gets executed in the next
// available _hot_ thread. after a thread finishes and no tasks are pending, it
// joins.
// TODO: take flag to keep threads hot after finishing a task to reduce overhead
// of launching a new thread with an optional timeout.
struct pool {
  using Task = std::function<void(std::uint16_t)>;

  pool(const pool &) = delete;
  pool(pool &&) = delete;

  explicit pool(std::uint16_t pool_size) {
    Assert(pool_size > 0, "pool size cannot be 0");
    Assert(pool_size <= thread::num_available_threads, "pool size", pool_size,
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

  std::size_t size() const { return _pool.size(); }

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
    logger::warn("cold start of worker", task_id);
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

// FIXME: turn logger into a class

// TODO: move to telemetry log level
struct scope_telemetry {
  explicit inline scope_telemetry(const char *msg)
      : callback([msg, start = std::chrono::high_resolution_clock::now()] {
          auto end = std::chrono::high_resolution_clock::now();
          logger::warn(msg, logger::duration(start, end));
        }) {}

  template <typename SizeType>
  inline scope_telemetry(const char *msg, SizeType size)
      : callback(
            [msg, size, start = std::chrono::high_resolution_clock::now()] {
              auto end = std::chrono::high_resolution_clock::now();
              logger::warn(msg, logger::throughput(start, end, size));
            }) {}
  // TODO add producer and consumer for progress
private:
  scope_dtor callback;
};

// FIXME: enable this macro with with telemetry flag
// #ifndef NDEBUG
// #define make_scope_timer(x) scope_timer x(#x)
// #else
// #define make_scope_timer(x) \
//   do { \ } while (0)
// #endif

template <typename Iterator> using range_pair = std::pair<Iterator, Iterator>;

#endif // CORE_HPP
