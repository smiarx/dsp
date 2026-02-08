#pragma once

#ifdef __clang__
// error with clang and -ffast-math
// only for clang > 1.18
#pragma clang diagnostic push
#if __clang__ > 1 || (__clang__ == 1 && __clang_major__ >= 18)
#pragma clang diagnostic ignored "-Wnan-infinity-disabled"
#endif

#include <gcem.hpp>

#ifdef __clang__
#pragma clang diagnostic pop
#endif
