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
#include <condition_variable>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <ostream>
#include <queue>
#include <sstream>
#include <string>

#include <thread>

namespace {
std::atomic_ulong thread_id_count = 0;
thread_local bool thread_id_assigned = false;
thread_local unsigned long thread_id = 0;
} // namespace

inline unsigned long get_thread_id() {
  if (not thread_id_assigned) {
    thread_id = thread_id_count.fetch_add(1, std::memory_order_relaxed);
    thread_id_assigned = true;
  }
  return thread_id;
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
static inline std::mutex mtx;
static inline std::queue<std::string> log_queue;
static inline std::mutex queue_mtx;
static inline std::condition_variable cv;
static inline std::thread logger_thread;
static inline bool stop_logger = false;

enum flags { none = 0x0, bold = 0b0001, ansi_colours = 0b0010 };

enum class LogLevel : unsigned {
  Off = 0,
  Error,
  Warn,
  Info,
  Debug

#ifndef LOG_LEVEL
#define LOG_LEVEL Info
#endif
};

static constexpr LogLevel current_level = LogLevel::LOG_LEVEL;

inline constexpr const char *get_level_str(LogLevel level) {
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

[[nodiscard]] static inline auto setup_logger() {
  std::cout << std::fixed << std::setprecision(3) << std::left;
  std::cerr << std::fixed << std::setprecision(3) << std::left;
  logger_thread = std::thread([] {
    while (true) {
      std::unique_lock<std::mutex> lock(queue_mtx);
      cv.wait(lock, [] { return !log_queue.empty() || stop_logger; });

      if (stop_logger && log_queue.empty()) {
        break;
      }

      std::string msg = std::move(log_queue.front());
      log_queue.pop();
      lock.unlock();
      std::cerr << msg.c_str() << std::endl;
    }
  });
  return scope_dtor([]() {
    {
      std::unique_lock<std::mutex> lock(queue_mtx);
      stop_logger = true;
    }
    cv.notify_all(); // Notify the logger thread to wake up and exit
    if (logger_thread.joinable()) {
      logger_thread.join();
    }
  });
}

template <typename Chrono> inline double duration(Chrono begin, Chrono end) {
  return std::chrono::duration<double>(end - begin).count();
}

template <typename Chrono>
inline std::string time_diff(Chrono begin, Chrono end,
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

inline std::string progress(std::size_t index, std::size_t size) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2) << std::left
      << ((double(index + 1) * 100 / size)) << "%";
  return oss.str();
}

template <LogLevel Level, typename... Args>
static inline void _log_impl(const char *color_code, Args &&...args) {
  if constexpr (Level <= current_level) {
    auto now = std::chrono::system_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                  now.time_since_epoch()) %
              1000;
    std::time_t tt = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf;
    localtime_r(&tt, &tm_buf);

    std::ostringstream oss;
    oss << color_code;
    oss << std::put_time(&tm_buf, "%Y-%m-%d %H:%M:%S") << '.'
        << std::setfill('0') << std::setw(3) << ms.count() << " ["
        << std::setfill('0') << std::setw(3) << get_thread_id() << "] "
        << get_level_str(Level) << " ";
    ((oss << std::forward<Args>(args) << " "), ...);
    oss << COL_RESET;
    {
      std::unique_lock<std::mutex> lock(queue_mtx);
      log_queue.emplace(oss.str());
    }
    cv.notify_one();
  }
}

template <typename... Args> inline void debug(Args &&...args) {
  _log_impl<LogLevel::Debug>(COL_CYAN, std::forward<Args>(args)...);
}

template <typename... Args> inline void info(Args &&...args) {
  _log_impl<LogLevel::Info>(COL_RESET, std::forward<Args>(args)...);
}

template <typename... Args> inline void warn(Args &&...args) {
  _log_impl<LogLevel::Warn>(COL_YELLOW, std::forward<Args>(args)...);
}

template <typename... Args> inline void error(Args &&...args) {
  _log_impl<LogLevel::Error>(COL_RED, std::forward<Args>(args)...);
  std::exit(EXIT_FAILURE);
}

template <typename... Args> inline void print(Args &&...args) {
  std::scoped_lock print_lock(mtx);
  ((std::cout << std::forward<Args>(args) << " "), ...);
  std::cout << COL_RESET << std::endl;
}

template <typename... Args> inline void printv(Args &&...args) {
  std::scoped_lock print_lock(mtx);
  ((std::cout << std::forward<Args>(args)), ...); // verbose
  std::cout << COL_RESET << std::endl;
}

template <LogLevel level = LogLevel::Info, typename... Args>
inline void print_thread(std::uint32_t thread_id, const char *color_code,
                         Args &&...args) {
  if constexpr (level <= current_level) {
    std::scoped_lock print_lock(mtx);
    std::cout << std::string((thread_id + 1), ' ') << color_code << std::setw(2)
              << thread_id << COL_RESET << " ";
    ((std::cout << std::forward<Args>(args) << " "), ...);
    std::cout << COL_RESET << std::endl;
  }
}

const auto logger_destoy = logger::setup_logger();
} // namespace logger

#endif // CORE_HPP
