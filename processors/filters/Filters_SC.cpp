#include "dsp/Context.h"
#include "dsp/VAFilters.h"
#include <SC_InterfaceTable.h>
#include <SC_PlugIn.h> // NOLINT
#include <SC_Rate.h>
#include <SC_Unit.h>

extern InterfaceTable *ft;

template <class Base, bool hasRes = false> struct Filter : public Unit {
    Base filter{};
    typename Base::State state{};
    float freq{0};
};

template <class Base> struct Filter<Base, true> : public Unit {
    Base filter{};
    typename Base::State state{};
    float freq{0};
    float res{0};
};

template <class Base, bool hasRes>
static void filterCtor(Filter<Base, hasRes> *unit);
template <class Base, bool hasRes>
static void filterNext(Filter<Base, hasRes> *unit, int inNumSamples);
template <class Base, bool hasRes>
static void filterNextA(Filter<Base, hasRes> *unit, int inNumSamples);

template <class Base, bool hasRes> void filterCtor(Filter<Base, hasRes> *unit)
{
    if (INRATE(1) == calc_FullRate) {
        auto &nextFunc = filterNextA<Base, hasRes>;
        SETCALC(nextFunc);
    } else {
        auto &nextFunc = filterNext<Base, hasRes>;
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
    unit->state = {};

    OUT0(0) = 0.f;
}

template <class Base, bool hasRes>
void filterNextA(Filter<Base, hasRes> *unit, int inNumSamples)
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

        float x = *in;
        unit->filter.process(dsp::Context(&x), unit->state);
        *out = x;
        ++in;
        ++out;
    }
}

template <class Base, bool hasRes>
void filterNext(Filter<Base, hasRes> *unit, int inNumSamples)
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
        float x = *in;
        unit->filter.process(dsp::Context(&x), unit->state);
        *out = x;
        ++in;
        ++out;
    }
}

#define CMA ,
#define DefineFilterUnit(name, Base, hasRes)                \
    (*ft->fDefineUnit)(#name, sizeof(Filter<Base, hasRes>), \
                       (UnitCtorFunc) & filterCtor<Base, hasRes>, 0, 0);
extern void loadFilters();
void loadFilters()
{
    DefineFilterUnit(OnePoleLP, dsp::va::OnePole<float CMA dsp::va::kLowPass>,
                     false);
    DefineFilterUnit(OnePoleHP, dsp::va::OnePole<float CMA dsp::va::kHighPass>,
                     false);
    DefineFilterUnit(OnePoleAP, dsp::va::OnePole<float CMA dsp::va::kAllPass>,
                     false);
    DefineFilterUnit(SVFLP, dsp::va::SVF<float CMA dsp::va::kLowPass>, true);
    DefineFilterUnit(SVFHP, dsp::va::SVF<float CMA dsp::va::kHighPass>, true);
    DefineFilterUnit(SVFAP, dsp::va::SVF<float CMA dsp::va::kAllPass>, true);
    DefineFilterUnit(SVFBP, dsp::va::SVF<float CMA dsp::va::kBandPass>, true);
    DefineFilterUnit(SVFNotch, dsp::va::SVF<float CMA dsp::va::kNotch>, true);
    DefineFilterUnit(LadderLP, dsp::va::Ladder<float CMA dsp::va::kLowPass>,
                     true);
    DefineFilterUnit(LadderHP, dsp::va::Ladder<float CMA dsp::va::kHighPass>,
                     true);
    DefineFilterUnit(LadderAP, dsp::va::Ladder<float CMA dsp::va::kAllPass>,
                     true);
}
