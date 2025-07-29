#pragma once

#include "../SC.h"
#include "dsp/cpu/defines.h"

void loadSprings();
#if DSP_X86_DISPATCH
void loadSpringsAVX2();
#endif

template <class Arch> struct SpringsUnit : public Unit {
    // processors::Springs *springs;
    Arch *springs;
};

template <class Arch>
static void springsNext(SpringsUnit<Arch> *unit, int inNumSamples);

template <class Arch> void Springs_Ctor(SpringsUnit<Arch> *unit)
{
    SETCALC(springsNext<Arch>);
    auto springs = (Arch *)RTAlloc(unit->mWorld, sizeof(Arch));
    new (springs) Arch();
    springs->prepare(SAMPLERATE, BUFLENGTH, [&unit](void *ptr, size_t len) {
        return RTRealloc(unit->mWorld, ptr, len);
    });
    unit->springs = springs;

    ClearUnitIfMemFailed(unit->springs);

    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
}

template <class Arch> void Springs_Dtor(SpringsUnit<Arch> *unit)
{
    unit->springs->free([&unit](void *ptr) { RTFree(unit->mWorld, ptr); });
    RTFree(unit->mWorld, unit->springs);
}

template <class Arch>
static void springsNext(SpringsUnit<Arch> *unit, int inNumSamples)
{
    float *in[2]  = {IN(9), IN(10)};
    float *out[2] = {OUT(0), OUT(1)};

    unit->springs->update(IN0(0), IN0(1), IN0(2), IN0(3), IN0(4), IN0(5),
                          IN0(6), IN0(7), IN0(8), inNumSamples);
    unit->springs->process(in, out, inNumSamples);
}

// load funcs
#define DefineSpringsUnit(Arch)                              \
    (*ft->fDefineUnit)("Springs", sizeof(SpringsUnit<Arch>), \
                       (UnitCtorFunc) & Springs_Ctor<Arch>,  \
                       (UnitDtorFunc) & Springs_Dtor<Arch>, 0);
