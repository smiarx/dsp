#include "dsp/cpu/defines.h"

#if DSP_X86_DISPATCH

#include "../Springs.h"
#include "../Springs_SC.h"

void loadSpringsAVX2() { DefineSpringsUnit(processors::Springs); }

#endif
