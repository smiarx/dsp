#include "dsp/cpu/infos.h"
#include <SC_InterfaceTable.h>

#include "adaptive/Adaptive_SC.h"
#include "springs/Springs_SC.h"
#include "tapedelay/TapeDelay_SC.h"
#include "voicebox/VoiceBox_SC.h"

extern void loadFilters();

InterfaceTable *ft;
PluginLoad(Processors)
{
    ft = inTable;

#ifdef DSP_X86_DISPATCH
    auto infos = dsp::cpu::getInfos();

    if (infos.avx2 && infos.fma3_sse42) {
        loadTapeDelayAVX2();
        loadSpringsAVX2();
        loadVoiceBox();
        loadFilters();
        // loadAdaptive();
        return;
    }
#endif

    loadTapeDelay();
    loadSprings();
    loadVoiceBox();
    loadFilters();
    // loadAdaptive();
}
