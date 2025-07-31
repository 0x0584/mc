
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

#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

#include <queue>
#include <vector>

#include <gperftools/profiler.h>

#include "thread.hpp"

using namespace std::chrono;

struct Foo {
  using duration_ms = std::chrono::duration<std::uint64_t, std::ratio<1, 1000>>;

  explicit Foo(duration_ms d) : sleep_ms(d) {
    logger::info("Foo(d)");
    id = task_id++;
  }

  Foo() { logger::info("Foo()"); }

  ~Foo() { logger::debug("~Foo()", id); }

  Foo(const Foo &other) : id(other.id), sleep_ms(other.sleep_ms) {
    logger::debug("Foo(const Foo &)", id);
  }

  Foo(Foo &&other)
      : id(std::move(other.id)), sleep_ms(std::move(other.sleep_ms)) {
    logger::debug("Foo(Foo &&):", id);
  }

  Foo &operator=(const Foo &rhs) {
    logger::debug("Foo::operator=(const Foo &)");
    sleep_ms = rhs.sleep_ms;
    id = rhs.id;
    return *this;
  }

  Foo &operator=(Foo &&rhs) {
    logger::debug(BG_YELLOW, COL_RED, "Foo::operator=(Foo &&)");
    sleep_ms = std::move(rhs.sleep_ms);
    id = std::move(rhs.id);
    return *this;
  }

  void operator()() {
    logger::info("Foo", id, "began executing..");
    std::this_thread::sleep_for(sleep_ms);
    logger::info("Foo", id, "finished.");
  }

private:
  static inline std::atomic_uint64_t task_id = 1;

  std::uint64_t id;
  duration_ms sleep_ms;
};

int main() {

  thread::pool pool(4);
  for (int i = 0; i < 10; ++i) {
    pool.exec(Foo(2 * 1000ms));
    pool.exec(Foo(3 * 300ms));
    pool.exec(Foo(5 * 800ms));
    pool.exec(Foo(9 * 100ms));
  }
}
