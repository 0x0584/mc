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

namespace memory {
std::pmr::unsynchronized_pool_resource *pool() {
  thread_local gc mem(__get_thread_id());
  return mem.get_pool();
}
} // namespace memory

scope_dtor logger::setup_logger() {
  {
    std::ostringstream oss;
    oss << "Logging started at " << std::chrono::system_clock::now() << '\n'
        << "Logging level is set to " << get_level_str(current_level) << '\n';
    locked_puts(oss.str(), stderr);
  }
  logger_thread = std::thread([] {
    std::pmr::list<log_entry> log_processing;
    while (true) {
      {
        std::unique_lock lk(queue_mtx);
        logger_cv.wait_for(lk, std::chrono::milliseconds(100), [] {
          bool ready = stop_logger || not log_queue.empty();
          if (std::scoped_lock flush_lock(flush_mtx); flush_logs && !ready) {
            flush_logs = false;
            flush_cv.notify_one();
          }
          return ready;
        });

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
      locked_puts(oss.str(), stderr);
      {
        std::scoped_lock lk(flush_mtx);
        if (flush_logs && log_processing.empty()) {
          flush_logs = false;
          flush_cv.notify_one();
        }
      }
    }
    std::ostringstream oss;
    oss << "Logging finished at " << std::chrono::system_clock::now() << '\n';
    locked_puts(oss.str(), stderr);
  });

  return scope_dtor([] {
    {
      std::scoped_lock lk(queue_mtx);
      stop_logger = true;
      logger_cv.notify_one();
    }
    if (logger_thread.joinable()) {
      logger_thread.join();
    }
  });
}
