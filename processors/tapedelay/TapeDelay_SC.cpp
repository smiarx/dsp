#include <SC_InterfaceTable.h>
#include <SC_PlugIn.h> // NOLINT
#include <SC_Unit.h>
#include <Unroll.h>
#include <cstddef>

#include "TapeDelay.h"

extern InterfaceTable *ft;

struct TapeDelay : public Unit {
    processors::TapeDelay *tapedelay;
};

static void TapeDelay_Ctor(TapeDelay *unit);
static void TapeDelay_Dtor(TapeDelay *unit);
static void tapeDelayNext(TapeDelay *unit, int inNumSamples);

void TapeDelay_Ctor(TapeDelay *unit)
{
    SETCALC(tapeDelayNext);

    unit->tapedelay = (processors::TapeDelay *)RTAlloc(
        unit->mWorld, sizeof(processors::TapeDelay));
    ClearUnitIfMemFailed(unit->tapedelay);
    new (unit->tapedelay) processors::TapeDelay();

    unit->tapedelay->prepare(SAMPLERATE, BUFLENGTH,
                             [&unit](void *ptr, size_t len) {
                                 return RTRealloc(unit->mWorld, ptr, len);
                             });

    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
}

void TapeDelay_Dtor(TapeDelay *unit)
{

    unit->tapedelay->free([&unit](void *ptr) { RTFree(unit->mWorld, ptr); });
    RTFree(unit->mWorld, unit->tapedelay);
}

void tapeDelayNext(TapeDelay *unit, int inNumSamples)
{
    float *in[2]  = {IN(8), IN(9)};
    float *out[2] = {OUT(0), OUT(1)};

    unit->tapedelay->update(IN0(0), IN0(1), IN0(2), IN0(3), IN0(4), IN0(5),
                            static_cast<processors::TapeDelay::Mode>(IN0(6)),
                            IN0(7), inNumSamples);
    unit->tapedelay->process(in, out, inNumSamples);
}

extern void loadTapeDelay();
void loadTapeDelay() { DefineDtorUnit(TapeDelay); }
