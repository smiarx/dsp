#include "dsp/VAFilters.h"
#include "dsp/Context.h"
#include "dsp/FastMath.h"
#include "dsp/Signal.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>

static constexpr float kSr      = 48000.f;
static constexpr float kNyquist = kSr / 2.f;

template <size_t N>
static void generatesine(dsp::fData<N> f, dsp::fData<N> phase,
                         dsp::fSample<N> *x, size_t length)
{
    // generate signal
    for (size_t n = 0; n < length; ++n) {
        for (size_t i = 0; i < N; ++i) {
            x[n][i] = sinf(dsp::constants<float>::pi * f[i] * n + phase[i]);
        }
    }
}

template <class Filter, size_t N>
static void filterArray(dsp::fData<N> cut, dsp::fSample<N> *x, size_t length)
{
    Filter filter;
    typename decltype(filter)::State filterstate;

    filter.setFreq(cut);
    filter.processBlock(dsp::Context(x, length), filterstate);
}

template <class Filter, size_t N>
static void filterResArray(dsp::fData<N> cut, dsp::fData<N> res,
                           dsp::fSample<N> *x, size_t length)
{
    Filter filter;
    typename decltype(filter)::State filterstate;

    filter.setFreq(cut, res);
    filter.processBlock(dsp::Context(x, length), filterstate);
}

template <size_t N>
static dsp::fSample<N> rms(dsp::fSample<N> *x, size_t length)
{
    dsp::fSample<N> rms = {0.f};
    for (size_t n = 0; n < length; ++n) {
        for (size_t i = 0; i < N; ++i) {
            rms[i] += x[n][i] * x[n][i];
        }
    }
    for (size_t i = 0; i < N; ++i) {
        rms[i] = sqrtf(rms[i] / length);
    }

    return rms;
}

TEST_CASE("VA Filters", "[dsp][va][onepole][lowpass]")
{

    // test filter reaction to DC signal
    SECTION("DC")
    {
        constexpr size_t N      = 1;
        constexpr size_t Length = 1 << 11;
        dsp::fSample<N> x[Length];

        auto cut = GENERATE(500.f / kNyquist, 1000.f / kNyquist,
                            2000.f / kNyquist, 15000.f / kNyquist);

#define TEST_DC(FILTER, TYPE, VAL)                                     \
    generatesine({0.f}, {dsp::constants<float>::pi / 2.f}, x, Length); \
    filterArray<FILTER<N, TYPE>>({cut}, x, Length);                    \
    for (size_t n = 6.f / cut; n < Length; ++n) {                      \
        REQUIRE_THAT(x[n][0], Catch::Matchers::WithinAbs(VAL, 1e-3));  \
    }

        TEST_DC(dsp::va::OnePole, dsp::va::kLowPass, 1.f);
        TEST_DC(dsp::va::SVF, dsp::va::kLowPass, 1.f);
        TEST_DC(dsp::va::Ladder, dsp::va::kLowPass, 1.f);

        TEST_DC(dsp::va::OnePole, dsp::va::kHighPass, 0.f);
        TEST_DC(dsp::va::SVF, dsp::va::kHighPass, 0.f);
        TEST_DC(dsp::va::Ladder, dsp::va::kHighPass, 0.f);
    }

    // test filter reaction to Cuttof frequency signal
    SECTION("Cutoff Frequency")
    {
        constexpr size_t N      = 1;
        constexpr size_t Length = 1 << 15;
        float cutGainm3db       = sqrtf(1.f / 2.f);
        float cutGainm12db      = 1 / 4.f;
        float power;
        dsp::fSample<N> x[Length];

        auto cut = GENERATE(take(6, random(0.f, 1.f)));

#define TEST_CUTOFF(FILTER, TYPE, GAIN)             \
    generatesine({cut}, {0.f}, x, Length);          \
    power = rms(x, Length)[0];                      \
    filterArray<FILTER<N, TYPE>>({cut}, x, Length); \
    REQUIRE_THAT(rms(x, Length)[0] / power,         \
                 Catch::Matchers::WithinAbs(GAIN, 5e-2));

        TEST_CUTOFF(dsp::va::OnePole, dsp::va::kLowPass, cutGainm3db);
        TEST_CUTOFF(dsp::va::SVF, dsp::va::kLowPass, cutGainm3db);
        TEST_CUTOFF(dsp::va::Ladder, dsp::va::kLowPass, cutGainm12db);

        TEST_CUTOFF(dsp::va::OnePole, dsp::va::kHighPass, cutGainm3db);
        TEST_CUTOFF(dsp::va::SVF, dsp::va::kHighPass, cutGainm3db);
        TEST_CUTOFF(dsp::va::Ladder, dsp::va::kHighPass, cutGainm12db);
    }

    // test filter reaction to high frequency signal
    SECTION("High Frequency")
    {
        constexpr size_t N         = 4;
        constexpr size_t Length    = 1 << 10;
        constexpr dsp::fData<N> fs = {0.97f, 0.964f, 0.9134f, 0.9273f};
        dsp::fSample<N> x[Length];
        dsp::fSample<N> powerOrig;
        dsp::fSample<N> power;

        auto cut = GENERATE(50.3f / kNyquist, 104.38f / kNyquist,
                            200.38f / kNyquist, 332.2f / kNyquist);

#define TEST_HIGHFREQ(FILTER, TYPE, VAL)                               \
    generatesine(fs, {0.f}, x, Length);                                \
    powerOrig = rms(x, Length);                                        \
    filterArray<FILTER<N, TYPE>>({cut, cut, cut, cut}, x, Length);     \
    power = rms(x, Length);                                            \
    for (size_t i = 0; i < N; ++i) {                                   \
        REQUIRE_THAT(power[i], Catch::Matchers::WithinAbs(VAL, 1e-2)); \
    }

        TEST_HIGHFREQ(dsp::va::OnePole, dsp::va::kLowPass, 0.f);
        TEST_HIGHFREQ(dsp::va::SVF, dsp::va::kLowPass, 0.f);
        TEST_HIGHFREQ(dsp::va::Ladder, dsp::va::kLowPass, 0.f);

        TEST_HIGHFREQ(dsp::va::OnePole, dsp::va::kHighPass, powerOrig[i]);
        TEST_HIGHFREQ(dsp::va::SVF, dsp::va::kHighPass, powerOrig[i]);
        TEST_HIGHFREQ(dsp::va::Ladder, dsp::va::kHighPass, powerOrig[i]);
    }

    // test filter resonance
    SECTION("Resonance")
    {
        constexpr size_t N      = 1;
        constexpr size_t Length = 1 << 15;
        float power;
        dsp::fSample<N> x[Length];

        auto cut = GENERATE(take(4, random(0.f, 0.7f)));
        auto R   = GENERATE(0.2f, take(4, random(0.08f, 2.f)));

        float cutGainm3db = 1.f / (2.f * R);

#define TEST_RES(FILTER, TYPE, GAIN)                        \
    generatesine({cut}, {0.f}, x, Length);                  \
    power = rms(x, Length)[0];                              \
    filterResArray<FILTER<N, TYPE>>({cut}, {R}, x, Length); \
    REQUIRE_THAT(rms(x, Length)[0] / power,                 \
                 Catch::Matchers::WithinAbs(GAIN, 5e-2));

        TEST_RES(dsp::va::SVF, dsp::va::kLowPass, cutGainm3db);
        TEST_RES(dsp::va::SVF, dsp::va::kHighPass, cutGainm3db);
        TEST_RES(dsp::va::SVF, dsp::va::kBandPass, 1.f);
    }
}
