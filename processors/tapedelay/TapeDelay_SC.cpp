#include "TapeDelay_SC.h"
#include "../SC.h"
#include "TapeDelay.h"
#include "dsp/cpu/defines.h"

template <class Arch> struct TapeDelayUnit : public Unit {
    Arch *tapedelay;
};

template <class Arch>
static void tapeDelayNext(TapeDelayUnit<Arch> *unit, int inNumSamples);

template <class Arch> void TapeDelay_Ctor(TapeDelayUnit<Arch> *unit)
{
    SETCALC(tapeDelayNext<Arch>);

    unit->tapedelay = (Arch *)RTAlloc(unit->mWorld, sizeof(Arch));
    ClearUnitIfMemFailed(unit->tapedelay);
    new (unit->tapedelay) Arch();

    unit->tapedelay->prepare(SAMPLERATE, BUFLENGTH,
                             [&unit](void *ptr, size_t len) {
                                 return RTRealloc(unit->mWorld, ptr, len);
                             });

    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
}

template <class Arch> void TapeDelay_Dtor(TapeDelayUnit<Arch> *unit)
{

    unit->tapedelay->free([&unit](void *ptr) { RTFree(unit->mWorld, ptr); });
    RTFree(unit->mWorld, unit->tapedelay);
}

template <class Arch>
void tapeDelayNext(TapeDelayUnit<Arch> *unit, int inNumSamples)
{
    float *in[2]  = {IN(8), IN(9)};
    float *out[2] = {OUT(0), OUT(1)};

    unit->tapedelay->update(IN0(0), IN0(1), IN0(2), IN0(3), IN0(4), IN0(5),
                            static_cast<typename Arch::Mode>(IN0(6)), IN0(7),
                            inNumSamples);
    unit->tapedelay->process(in, out, inNumSamples);
}

// load funcs
#define DefineTapeDelayUnit(Arch)                                \
    (*ft->fDefineUnit)("TapeDelay", sizeof(TapeDelayUnit<Arch>), \
                       (UnitCtorFunc) & TapeDelay_Ctor<Arch>,    \
                       (UnitDtorFunc) & TapeDelay_Dtor<Arch>, 0);

#if DSP_COMPILE_AVX2
#define LOADFUNC loadTapeDelayAVX2
#else
#define LOADFUNC loadTapeDelay
#endif
void LOADFUNC() { DefineTapeDelayUnit(processors::TapeDelay); }
