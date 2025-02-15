#include "../../Context.h"
#include "../../VAFilters.h"
#include "SC_PlugIn.h"
#include "SC_Unit.h"

static InterfaceTable *ft;

template <class Base, bool hasRes = false> struct Filter : public Unit {
    Base filter;
    typename Base::State state;
    float freq{0};
};

template <class Base> struct Filter<Base, true> : public Unit {
    Base filter;
    typename Base::State state;
    float freq{0};
    float res{0};
};

template <class Base, bool hasRes> void Filter_Ctor(Filter<Base, hasRes> *unit);
template <class Base, bool hasRes>
void Filter_next(Filter<Base, hasRes> *unit, int inNumSamples);
template <class Base, bool hasRes>
void Filter_next_a(Filter<Base, hasRes> *unit, int inNumSamples);

template <class Base, bool hasRes> void Filter_Ctor(Filter<Base, hasRes> *unit)
{
    if (INRATE(1) == calc_FullRate) {
        auto &nextFunc = Filter_next_a<Base, hasRes>;
        SETCALC(nextFunc);
    } else {
        auto &nextFunc = Filter_next<Base, hasRes>;
        SETCALC(nextFunc);
    }

    unit->freq = IN0(1);
    float freq = unit->freq * SAMPLEDUR * 2.f;
    if constexpr (hasRes) {
        unit->res    = IN0(2);
        unit->filter = Base({freq}, {unit->res});
    } else {
        unit->filter = Base({freq});
    }
    unit->state = {{{0.f}}};

    OUT0(0) = 0.f;
}

template <class Base, bool hasRes>
void Filter_next_a(Filter<Base, hasRes> *unit, int inNumSamples)
{
    if (unit->freq != IN0(1)) {
        unit->freq = IN0(1);
    } else {
        if constexpr (hasRes) {
            float res = IN0(2);
            if (unit->res != res) {
                unit->filter.setRes({res});
                unit->res = res;
            }
        }
    }

    auto *in  = IN(0);
    auto *out = OUT(0);
    for (int n = 0; n < inNumSamples; ++n) {
        float freq = IN(1)[n] * SAMPLEDUR * 2.f;
        if constexpr (hasRes) {
            float res    = IN0(2);
            unit->filter = Base({freq}, {res});
        } else {
            unit->filter = Base({freq});
        }

        dsp::Signal<1> x = {*in};
        unit->filter.process(dsp::Context(&x), unit->state);
        *out = x[0];
        ++in;
        ++out;
    }
}

template <class Base, bool hasRes>
void Filter_next(Filter<Base, hasRes> *unit, int inNumSamples)
{
    if (unit->freq != IN0(1)) {
        float freq = IN0(1) * SAMPLEDUR * 2.f;
        if constexpr (hasRes) {
            float res    = IN0(2);
            unit->filter = Base({freq}, {res});
            unit->res    = res;
        } else {
            unit->filter = Base({freq});
        }
        unit->freq = IN0(1);
    } else {
        if constexpr (hasRes) {
            float res = IN0(2);
            if (unit->res != res) {
                unit->filter.setRes({res});
                unit->res = res;
            }
        }
    }

    auto *in  = IN(0);
    auto *out = OUT(0);
    for (int n = 0; n < inNumSamples; ++n) {
        dsp::Signal<1> x = {*in};
        unit->filter.process(dsp::Context(&x), unit->state);
        *out = x[0];
        ++in;
        ++out;
    }
}

#define CMA ,
#define DefineFilterUnit(name, Base, hasRes)                \
    (*ft->fDefineUnit)(#name, sizeof(Filter<Base, hasRes>), \
                       (UnitCtorFunc) & Filter_Ctor<Base, hasRes>, 0, 0);
PluginLoad(SCTapeDelay)
{
    ft = inTable; // store pointer to InterfaceTable
    DefineFilterUnit(OnePoleLP, dsp::OnePole<1 CMA dsp::LowPass>, false);
    DefineFilterUnit(OnePoleHP, dsp::OnePole<1 CMA dsp::HighPass>, false);
    DefineFilterUnit(OnePoleAP, dsp::OnePole<1 CMA dsp::AllPass>, false);
    DefineFilterUnit(SVFLP, dsp::SVF<1 CMA dsp::LowPass>, true);
    DefineFilterUnit(SVFHP, dsp::SVF<1 CMA dsp::HighPass>, true);
    DefineFilterUnit(SVFAP, dsp::SVF<1 CMA dsp::AllPass>, true);
    DefineFilterUnit(SVFBP, dsp::SVF<1 CMA dsp::BandPass>, true);
    DefineFilterUnit(SVFNotch, dsp::SVF<1 CMA dsp::Notch>, true);
    DefineFilterUnit(LadderLP, dsp::Ladder<1 CMA dsp::LowPass>, true);
    DefineFilterUnit(LadderHP, dsp::Ladder<1 CMA dsp::HighPass>, true);
    DefineFilterUnit(LadderAP, dsp::Ladder<1 CMA dsp::AllPass>, true);
}
