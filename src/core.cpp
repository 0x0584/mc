#include "core.hpp"

#ifndef THREADS_PER_CORE
#define THREADS_PER_CORE 2
#endif

namespace thread {
const std::uint16_t threads_per_core = THREADS_PER_CORE;
const std::uint16_t num_available_threads =
    std::thread::hardware_concurrency() * threads_per_core;
} // namespace thread

namespace {
std::atomic_uint16_t __thread_id_count = 0;
thread_local bool __thread_id_assigned = false;
thread_local std::uint16_t __thread_id = 0;
} // namespace

std::uint16_t __get_thread_id() {
  if (not __thread_id_assigned) {
    __thread_id = __thread_id_count.fetch_add(1, std::memory_order_relaxed);
    __thread_id_assigned = true;
  }
  return __thread_id;
}

scope_dtor logger::setup_logger() {
  std::cout << std::fixed << std::setprecision(3) << std::left;
  std::cerr << std::fixed << std::setprecision(3) << std::left;
  std::cerr << "Logging started at " << std::chrono::system_clock::now() << '\n'
            << "Logging level is set to " << get_level_str(current_level)
            << '\n';
  logger_thread = std::thread([] {
    std::pmr::list<log_entry> log_processing;

    while (true) {
      {
        std::unique_lock queue_lock(queue_mtx);
        logger_cv.wait_for(queue_lock, std::chrono::milliseconds(100),
                           [] { return stop_logger || not log_queue.empty(); });

        if (log_queue.empty()) {
          if (log_processing.empty() && stop_logger) {
            break;
          } else {
            continue;
          }
        }

        log_processing.splice(log_processing.begin(), log_queue);
      }

      // XXX: sort the logs by id

      log_processing.sort(std::less<log_entry>());

      // XXX: only print the current available logs
      // XXX: push the rest to the batch

      std::ostringstream oss;
      while (not log_processing.empty() &&
             log_processing.front().id() == print_log_id) {
        oss << std::move(log_processing.front().message());
        log_processing.pop_front();
        print_log_id++;
      }

      std::scoped_lock print_lock(print_mtx);
      std::cerr << std::move(oss.str());
    }
    std::cerr << "Logging finished at " << std::chrono::system_clock::now()
              << '\n';
  });

  return scope_dtor([] {
    {
      std::scoped_lock lock(queue_mtx);
      stop_logger = true;
      logger_cv.notify_one();
    }
    if (logger_thread.joinable()) {
      logger_thread.join();
    }
    std::cerr << std::flush;
  });
}
