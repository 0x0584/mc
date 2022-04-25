#ifndef UTILS_HPP
#define UTILS_HPP

#include <cassert>
#include <cstring>

#include <algorithm>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <execution>
#include <memory_resource>
#include <mutex>
#include <shared_mutex>
#include <thread>
#include <utility>

#include <iostream>
#include <sstream>
#include <stdexcept>

#define LITERAL(expr) #expr

#define TIMER(name) name##_timer
#define TIMER_NAME(name) LITERAL(name##_timer)

//#ifndef NDEBUG
#define MAKE_TIMER(name) timer TIMER(name)(TIMER_NAME(name))
#define STOP_TIMER(name) TIMER(name).stop();
//#else
// #define MAKE_TIME(name) void *name##timer
// #define STOP_TIMER(name)
// #endif

#ifndef THREADS_PER_CORE
#define THREADS_PER_CORE 8
#endif

namespace mc {
using size_type = unsigned;

extern constexpr std::uint16_t threads_per_core = THREADS_PER_CORE;
extern const std::uint16_t num_threads =
    std::thread::hardware_concurrency() * threads_per_core;

extern constexpr std::uint16_t default_thread = 0;

class timer {
  const char *name;
  std::chrono::time_point<std::chrono::high_resolution_clock> begin, end;
  bool set = false;

public:
  inline timer(const char *name = "timer")
      : name(name), begin(std::chrono::high_resolution_clock::now()) {}

  ~timer() {
    if (not set) {
      end = std::chrono::high_resolution_clock::now();
    }
    // clang-format off
    std::cerr << name  << " "
              << std::chrono::duration<double>(end - begin).count() << "s ("
              << std::chrono::duration<double, std::milli>(end - begin) .count() << "ms)\n";
    // clang-format on
  }

  inline void stop() {
    end = std::chrono::high_resolution_clock::now();
    set = true;
  }
};

template <typename Callable> struct end_of_scope_executor {
  end_of_scope_executor(const end_of_scope_executor &) = delete;
  end_of_scope_executor(end_of_scope_executor &&) = delete;

  inline explicit end_of_scope_executor(Callable &&fn)
      : callback(std::forward<Callable>(fn)) {}
  inline ~end_of_scope_executor() { callback(); }

  end_of_scope_executor &operator=(const end_of_scope_executor &) = delete;
  end_of_scope_executor &operator=(end_of_scope_executor &&) = delete;

private:
  Callable callback;
};
} // namespace mc
#endif // UTILS_HPP
