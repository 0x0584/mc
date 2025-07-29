// thread.hpp
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

#ifndef THREAD_HPP
#define THREAD_HPP

#ifndef THREADS_PER_CORE
#define THREADS_PER_CORE 2
#endif

#include <cassert>
#include <cstdint>
#include <queue>
#include <thread>

#include "core.hpp"

namespace thread {
const std::uint32_t threads_per_core = THREADS_PER_CORE;
const std::uint32_t num_threads =
    std::thread::hardware_concurrency() * threads_per_core;

template <typename Task> struct pool {
  static inline const std::uint16_t max_num_threads = 255;

  pool(const pool &) = delete;
  pool(pool &&) = delete;

  explicit pool(std::uint16_t pool_size = 1) {
    assert(pool_size > 0);
    assert(pool_size <= num_threads);
    _pool.resize(pool_size);
    for (std::uint16_t thread_id = 0u; thread_id < _pool.size(); ++thread_id) {
      _available.push(thread_id);
    }
  }

  ~pool() {
    logger::debug("~pool()");
    join();
  }

  pool &operator=(const pool &) = delete;
  pool &operator=(pool &&) = delete;

  void exec(Task &&task) {
    std::uint16_t thread_id;
    if (std::unique_lock<std::mutex> pool_lock(_pool_mtx); _available.empty()) {
      _pending_tasks.emplace(std::move(task));
      pool_lock.unlock();
      logger::debug("pool is full!");
      return;
    } else {
      thread_id = _available.front();
      _available.pop();
    }
    logger::debug(thread_id, "got task..");
    _pool[thread_id] = std::thread(
        [this, thread_id](Task callback) mutable {
          callback();
          _process_pending(thread_id);
          _available.emplace(thread_id);
          logger::info(thread_id, "is available");
        },
        std::forward<Task>(task));
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
      logger::debug(thread_id, "is processing a pending task");
      pending_task();
    }
  }

  std::queue<Task> _pending_tasks;
  std::queue<std::uint16_t> _available;
  std::vector<std::thread> _pool;
  std::mutex _pool_mtx;
};
} // namespace thread
#endif // THREAD_HPP
