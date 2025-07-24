// log.hpp
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

#ifndef LOG_HPP
#define LOG_HPP

#define COL_GREEN "\x1b[32m"
#define COL_MAGENTA "\x1b[35m"
#define COL_CYAN "\x1b[36m"
#define COL_RESET "\x1b[0m"
#define COL_BOLD "\x1b[1m"

#include <iomanip>
#include <iostream>
#include <mutex>
#include <ostream>
#include <sstream>
#include <string>

namespace mc {
struct log {
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

  template <typename... Args> static inline void info(Args &&...args) {
    std::scoped_lock print_lock(mtx);
    std::cerr << "> ";
    ((std::cerr << std::forward<Args>(args) << " "), ...);
    std::cerr << COL_RESET << std::endl;
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

  template <typename... Args> static inline void error(Args &&...args) {
    std::scoped_lock print_lock(mtx);
    ((std::cerr << std::forward<Args>(args) << " "), ...);
    std::cerr << COL_RESET << std::endl;
    std::exit(EXIT_FAILURE);
  }

  template <typename... Args>
  static inline void print_thread(std::uint32_t thread_id, Args &&...args) {
    std::scoped_lock print_lock(mtx);
    std::cout << std::string((thread_id + 1), ' ') << std::to_string(thread_id)
              << " ";
    ((std::cout << std::forward<Args>(args) << " "), ...);
    std::cout << COL_RESET << std::endl;
  }

private:
  static inline std::mutex mtx;
};
} // namespace mc
#endif // LOG_HPP
