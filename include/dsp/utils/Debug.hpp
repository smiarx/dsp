#pragma once

#if !defined(NDEBUG) && !defined(NDEBUGLOG)
#include <iostream>
#define DEBUGLOG std::cout
#else
#define DEBUGLOG \
    if constexpr (false) std::cout
#endif
