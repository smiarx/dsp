#pragma once

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#if __clang__ > 1 || (__clang__ == 1 && __clang_major__ >= 18)
// error with clang and -ffast-math
// only for clang > 1.18
#pragma clang diagnostic ignored "-Wnan-infinity-disabled"
#endif
#endif

#include <gcem.hpp>

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
