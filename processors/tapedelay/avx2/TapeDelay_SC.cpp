#include "dsp/cpu/defines.h"

#if DSP_X86_DISPATCH

#include "../TapeDelay.h"
#include "../TapeDelay_SC.h"

void loadTapeDelayAVX2() { DefineTapeDelayUnit(processors::TapeDelay); }
#endif
