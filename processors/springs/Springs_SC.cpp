#include "SC_PlugIn.h"
#include "Springs.h"

extern InterfaceTable *ft;

struct Springs : public Unit {
    processors::Springs *springs;
};

void Springs_Ctor(Springs *unit);
void Springs_Dtor(Springs *unit);
void Springs_next(Springs *unit, int inNumSamples);

void Springs_Ctor(Springs *unit)
{
    SETCALC(Springs_next);

    unit->springs = (processors::Springs *)RTAlloc(unit->mWorld,
                                                   sizeof(processors::Springs));
    ClearUnitIfMemFailed(unit->springs);
    new (unit->springs) processors::Springs();

    unit->springs->prepare(SAMPLERATE, BUFLENGTH,
                           [&unit](void *ptr, size_t len) {
                               return RTRealloc(unit->mWorld, ptr, len);
                           });

    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
}

void Springs_Dtor(Springs *unit)
{
    unit->springs->free([&unit](void *ptr) { RTFree(unit->mWorld, ptr); });
    RTFree(unit->mWorld, unit->springs);
}

void Springs_next(Springs *unit, int inNumSamples)
{
    float *in[2]  = {IN(9), IN(10)};
    float *out[2] = {OUT(0), OUT(1)};

    unit->springs->update(IN0(0), IN0(1), IN0(2), IN0(3), IN0(4), IN0(5),
                          IN0(6), IN0(7), IN0(8), inNumSamples);
    unit->springs->process(in, out, inNumSamples);
}

void LoadSprings() { DefineDtorUnit(Springs); }
