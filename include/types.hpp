#ifndef TYPES_HPP
#define TYPES_HPP

#include <stdint.h>

namespace mc {

namespace type {
typedef uint32_t key;
typedef uint64_t offset;
typedef key degree;
typedef uint16_t colour;
typedef uint32_t timestamp;
} // namespace type

using namespace type;
} // namespace mc
#endif
