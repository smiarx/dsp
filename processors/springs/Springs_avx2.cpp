#include "Springs_avx2.h"

#if DSP_X86_DISPATCH
#pragma GCC target "avx2"
#pragma GCC target "fma"

#define SPRINGS_NAMESPACE avx2
#include "Springs.cpp"
#endif
