#include "VoiceBox_SC.h"
#include "VoiceBox.h"

struct VoiceBoxUnit : public Unit {
    processors::VoiceBox *voicebox_;
};

static void voiceBoxNext(VoiceBoxUnit *unit, int inNumSamples);

static void voiceBoxCtor(VoiceBoxUnit *unit)
{
    SETCALC(voiceBoxNext);
    auto *voicebox = (processors::VoiceBox *)RTAlloc(
        unit->mWorld, sizeof(processors::VoiceBox));
    new (voicebox) processors::VoiceBox();
    unit->voicebox_ = voicebox;

    ClearUnitIfMemFailed(voicebox);
    ZOUT0(0) = 0.f;
}

static void voiceBoxDtor(VoiceBoxUnit *unit)
{
    RTFree(unit->mWorld, unit->voicebox_);
}

void voiceBoxNext(VoiceBoxUnit *unit, int inNumSamples)
{
    unit->voicebox_->update(IN0(1), IN0(2), IN0(3));
    unit->voicebox_->process(&IN(0), &OUT(0), inNumSamples);
}

void loadVoiceBox()
{
    (*ft->fDefineUnit)("VoiceBox", sizeof(VoiceBoxUnit),
                       (UnitCtorFunc)&voiceBoxCtor, (UnitDtorFunc)&voiceBoxDtor,
                       0);
}
