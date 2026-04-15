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
    unit->voicebox_->prepare(static_cast<float>(SAMPLERATE));

    ClearUnitIfMemFailed(voicebox);
    ZOUT0(0) = 0.f;
}

static void voiceBoxDtor(VoiceBoxUnit *unit)
{
    RTFree(unit->mWorld, unit->voicebox_);
}

void voiceBoxNext(VoiceBoxUnit *unit, int inNumSamples)
{
    unit->voicebox_->update(IN0(1), IN0(2));
    unit->voicebox_->process(&IN(0), &OUT(0), inNumSamples);
}

//===========================================================

#include <dsp/analysis/ZeroCrossing.hpp>
#include <dsp/filters/OnePole.hpp>

struct ZCUnit : public Unit {
    dsp::OnePole<float> op_;
    decltype(op_)::State opState_;
    dsp::ZeroCrossing<float> zc_;
};

static void zcNext(ZCUnit *unit, int inNumSamples)
{
    unit->op_.setRate(0.1f, 1.f / static_cast<float>(IN0(1) * SAMPLERATE));
    for (int n = 0; n < inNumSamples; ++n) {
        float v = unit->zc_.process(dsp::Context(&IN(0)[n])) ? 1.f : 0.f;
        unit->op_.process(dsp::Context(&v), unit->opState_);
        OUT(0)[n] = v;
    };
}

static void zcCtor(ZCUnit *unit)
{
    SETCALC(zcNext);
    unit->op_      = {};
    unit->opState_ = {};
    unit->zc_      = {};
}

//==================================================

#include "PitchTracker.hpp"

struct PitchTrackerUnit : public Unit {
    processors::PitchTracker *pitchTracker_;
};

static void pitchTrackerNext(PitchTrackerUnit *unit, int inNumSamples);

static void pitchTrackerCtor(PitchTrackerUnit *unit)
{
    SETCALC(pitchTrackerNext);
    auto *pitchTracker = (processors::PitchTracker *)RTAlloc(
        unit->mWorld, sizeof(processors::PitchTracker));
    new (pitchTracker) processors::PitchTracker();
    unit->pitchTracker_ = pitchTracker;
    unit->pitchTracker_->prepare(static_cast<float>(SAMPLERATE));

    ClearUnitIfMemFailed(pitchTracker);
    ZOUT0(0) = 0.f;
}

static void pitchTrackerDtor(PitchTrackerUnit *unit)
{
    RTFree(unit->mWorld, unit->pitchTracker_);
}

void pitchTrackerNext(PitchTrackerUnit *unit, int inNumSamples)
{
    unit->pitchTracker_->process(&IN(0), &OUT(0), inNumSamples);
    for (int i = 0; i < inNumSamples; ++i) {
        OUT(0)[i] *= static_cast<float>(SAMPLERATE / 2.f);
    }
}

//==================================================

void loadVoiceBox()
{
    (*ft->fDefineUnit)("VoiceBox", sizeof(VoiceBoxUnit),
                       (UnitCtorFunc)&voiceBoxCtor, (UnitDtorFunc)&voiceBoxDtor,
                       0);
    (*ft->fDefineUnit)("ZC", sizeof(ZCUnit), (UnitCtorFunc)&zcCtor, nullptr, 0);

    (*ft->fDefineUnit)("PitchTracker", sizeof(PitchTrackerUnit),
                       (UnitCtorFunc)&pitchTrackerCtor,
                       (UnitDtorFunc)&pitchTrackerDtor, 0);
}
