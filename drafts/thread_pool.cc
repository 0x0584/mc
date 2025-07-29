
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

#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

#include <cassert>
#include <cstdint>

#include <atomic>
#include <mutex>
#include <thread>

#include <queue>
#include <vector>

#include <gperftools/profiler.h>

namespace {
std::atomic_ulong thread_id_count = 0;
thread_local bool thread_id_assigned = false;
thread_local unsigned long thread_id = 0;
} // namespace

static inline unsigned long get_thread_id() {
  if (not thread_id_assigned) {
    thread_id = thread_id_count.fetch_add(1, std::memory_order_relaxed);
    thread_id_assigned = true;
  }
  return thread_id;
}

struct log {

  enum class LogLevel : unsigned { Off = 0, Error, Warn, Info, Debug };

  static constexpr LogLevel current_level = LogLevel::Debug;

  static constexpr const char *get_level_str(LogLevel level) {
    switch (level) {
    case LogLevel::Debug:
      return "DEBUG";
    case LogLevel::Info:
      return "INFO ";
    case LogLevel::Warn:
      return "WARN ";
    case LogLevel::Error:
      return "ERROR";
    default:
      throw std::logic_error("unknown LogLevel!");
    }
  }

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

  template <LogLevel Level> static constexpr std::string timestamp() {
    auto now = std::chrono::system_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                  now.time_since_epoch()) %
              1000;
    std::time_t tt = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf;
    localtime_r(&tt, &tm_buf);

    std::ostringstream oss;
    oss << std::put_time(&tm_buf, "%Y-%m-%d %H:%M:%S") << '.'
        << std::setfill('0') << std::setw(3) << ms.count() << " ["
        << std::setfill('0') << std::setw(3) << get_thread_id() << "] "
        << get_level_str(Level) << " ";

    return oss.str();
  }

  template <LogLevel Level, typename... Args>
  static inline void _log_impl(const char *color_code, Args &&...args) {
    if constexpr (Level <= current_level) {
      std::string stamp = timestamp<Level>();
      std::scoped_lock print_lock(mtx);
      std::cerr << color_code << stamp.c_str();
      ((std::cerr << std::forward<Args>(args) << " "), ...);
      std::cerr << COL_RESET << std::endl;
    }
  }

  template <typename... Args> static inline void debug(Args &&...args) {
    _log_impl<LogLevel::Debug>(COL_CYAN, std::forward<Args>(args)...);
  }

  template <typename... Args> static inline void info(Args &&...args) {
    _log_impl<LogLevel::Info>(COL_RESET, std::forward<Args>(args)...);
  }

  template <typename... Args> static inline void warn(Args &&...args) {
    _log_impl<LogLevel::Warn>(COL_YELLOW, std::forward<Args>(args)...);
  }

  template <typename... Args> static inline void error(Args &&...args) {
    _log_impl<LogLevel::Error>(COL_RED, std::forward<Args>(args)...);
    std::exit(EXIT_FAILURE);
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

  template <LogLevel level = LogLevel::Info, typename... Args>
  static inline void print_thread(std::uint32_t thread_id,
                                  const char *color_code, Args &&...args) {
    if constexpr (level <= current_level) {
      std::scoped_lock print_lock(mtx);
      std::cout << std::string((thread_id + 1), ' ') << color_code
                << std::setw(2) << thread_id << COL_RESET << " ";
      ((std::cout << std::forward<Args>(args) << " "), ...);
      std::cout << COL_RESET << std::endl;
    }
  }

private:
  static inline std::mutex mtx;
};

namespace thread {
const std::uint16_t num_threads = 8;

template <typename Task> struct pool {
  static inline const std::uint16_t max_num_threads = 255;

  pool(const pool &) = delete;
  pool(pool &&) = delete;

  explicit pool(std::uint16_t pool_size = 1) {
    assert(pool_size > 0);
    assert(pool_size <= std::min(max_num_threads, num_threads));
    _pool.resize(pool_size);
    for (std::uint16_t thread_id = 0u; thread_id < _pool.size(); ++thread_id) {
      _available.push(thread_id);
      log::info(thread_id, "is available");
    }
  }

  ~pool() {
    log::debug("~pool()");
    join();
  }

  pool &operator=(const pool &) = delete;
  pool &operator=(pool &&) = delete;

  void exec(Task &&task) {
    std::uint16_t thread_id;
    if (std::unique_lock<std::mutex> pool_lock(_pool_mtx); _available.empty()) {
      _pending_tasks.emplace(std::move(task));
      pool_lock.unlock();
      log::info("pool is full!");
      return;
    } else {
      thread_id = _available.front();
      _available.pop();
    }
    log::info(thread_id, "got task..");
    _pool[thread_id] = std::thread(
        [this, thread_id](Task callback) mutable {
          callback();
          _process_pending(thread_id);
          _available.emplace(thread_id);
          log::info(thread_id, "is available");
        },
        std::forward<Task>(task));
  }

  void join() {
    log::info("joining threads..");
    for (auto &task : _pool) {
      if (task.joinable()) {
        task.join();
      }
    }
    log::info("all threads joined.");
  }

private:
  void _process_pending(std::uint16_t thread_id) {
    while (true) {
      std::unique_lock<std::mutex> pool_lock(_pool_mtx);
      if (_pending_tasks.empty()) {
        break;
      }
      Task pending_task = std::move(_pending_tasks.front());
      _pending_tasks.pop();
      pool_lock.unlock();
      log::info(thread_id, "is processing a pending task");
      pending_task();
    }
  }

  std::queue<Task> _pending_tasks;
  std::queue<std::uint16_t> _available;
  std::vector<std::thread> _pool;
  std::mutex _pool_mtx;
};
} // namespace thread

using namespace std::chrono;

struct Foo {
  using duration_ms = std::chrono::duration<std::uint64_t, std::ratio<1, 1000>>;

  explicit Foo(duration_ms d) : sleep_ms(d) {
    log::info("Foo(d)");
    id = task_id++;
  }

  Foo() { log::info("Foo()"); }

  ~Foo() { log::debug("~Foo()", id); }

  Foo(const Foo &other) : id(other.id), sleep_ms(other.sleep_ms) {
    log::debug("Foo(const Foo &)", id);
  }

  Foo(Foo &&other)
      : id(std::move(other.id)), sleep_ms(std::move(other.sleep_ms)) {
    log::debug("Foo(Foo &&):", id);
  }

  Foo &operator=(const Foo &rhs) {
    log::debug("Foo::operator=(const Foo &)");
    sleep_ms = rhs.sleep_ms;
    id = rhs.id;
    return *this;
  }

  Foo &operator=(Foo &&rhs) {
    log::debug(BG_YELLOW, COL_RED, "Foo::operator=(Foo &&)");
    sleep_ms = std::move(rhs.sleep_ms);
    id = std::move(rhs.id);
    return *this;
  }

  void operator()() {
    log::info("Foo", id, "began executing..");
    std::this_thread::sleep_for(sleep_ms);
    log::info("Foo", id, "finished.");
  }

private:
  static inline std::atomic_uint64_t task_id = 1;

  std::uint64_t id;
  duration_ms sleep_ms;
};

int main() {
  log::setup_logger();

  thread::pool<Foo> pool(4);

  for (int i = 0; i < 1; ++i) {
    pool.exec(Foo(4 * 1000ms));
    pool.exec(Foo(3 * 1000ms));
    pool.exec(Foo(2 * 1000ms));
    pool.exec(Foo(1000ms));
  }
}
