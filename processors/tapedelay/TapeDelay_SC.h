#pragma once

#include "dsp/cpu/defines.h"

void loadTapeDelay();
#if DSP_X86_DISPATCH
void loadTapeDelayAVX2();
#endif
