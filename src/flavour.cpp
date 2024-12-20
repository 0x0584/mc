// flavour.cpp
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

#include "flavour.hpp"
namespace mc {
std::ostream &operator<<(std::ostream &os, flavour algo_type) {
  if (algo_type == flavour::exact) {
    return os << "Exact Algorithm";
  } else if (algo_type == flavour::heuristic) {
    return os << "Heuristic Algorithm";
  } else {
    return os << "Hybrid Algorithm";
  }
}
}