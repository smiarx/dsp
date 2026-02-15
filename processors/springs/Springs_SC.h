#pragma once

#include "dsp/cpu/defines.h"

void loadSprings();
#if DSP_X86_DISPATCH
void loadSpringsAVX2();
#endif
