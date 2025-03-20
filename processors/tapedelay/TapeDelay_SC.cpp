#include "SC_PlugIn.h"
#include "SC_Unit.h"

#include "TapeDelay.h"

extern InterfaceTable *ft;

struct TapeDelay : public Unit {
    processors::TapeDelay *tapedelay;
};

void TapeDelay_Ctor(TapeDelay *unit);
void TapeDelay_Dtor(TapeDelay *unit);
void TapeDelay_next(TapeDelay *unit, int inNumSamples);

void TapeDelay_Ctor(TapeDelay *unit)
{
    SETCALC(TapeDelay_next);

    unit->tapedelay = (processors::TapeDelay *)RTAlloc(
        unit->mWorld, sizeof(processors::TapeDelay));
    ClearUnitIfMemFailed(unit->tapedelay);
    new (unit->tapedelay) processors::TapeDelay();

    unit->tapedelay->setSampleRate(SAMPLERATE);

    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
}

void TapeDelay_Dtor(TapeDelay *unit) { RTFree(unit->mWorld, unit->tapedelay); }

void TapeDelay_next(TapeDelay *unit, int inNumSamples)
{
    float *in[2]  = {IN(8), IN(9)};
    float *out[2] = {OUT(0), OUT(1)};

    unit->tapedelay->update(IN0(0), IN0(1), IN0(2), IN0(3), IN0(4), IN0(5),
                            static_cast<processors::TapeDelay::Mode>(IN0(6)),
                            IN0(7), inNumSamples);
    unit->tapedelay->process(in, out, inNumSamples);
}

void LoadTapeDelay() { DefineDtorUnit(TapeDelay); }
