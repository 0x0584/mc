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
#define THREADS_PER_CORE 8
#endif

#include <cstdint>
#include <thread>

#include "log.hpp"

namespace thread {
const std::uint32_t threads_per_core = THREADS_PER_CORE;
const std::uint32_t num_threads =
    std::thread::hardware_concurrency() * threads_per_core;

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
}; // namespace thread
#endif // THREAD_HPP