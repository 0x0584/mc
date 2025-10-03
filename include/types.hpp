#ifndef TYPES_HPP
#define TYPES_HPP

#include <numeric>
#include <stdint.h>

namespace mc {

namespace type {
typedef uint32_t key;
typedef uint64_t offset;
typedef key degree;
typedef uint16_t colour;
typedef uint32_t timestamp;

static inline constexpr key key_npos = std::numeric_limits<key>::max();
static inline constexpr key offset_npos = std::numeric_limits<offset>::max();
} // namespace type

using namespace type;
} // namespace mc
#endif
