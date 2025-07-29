#include "dsp/VAFilters.h"
#include "dsp/Context.h"
#include "dsp/FastMath.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>

static constexpr float kSr      = 48000.f;
static constexpr float kNyquist = kSr / 2.f;

template <typename T>
static void generatesine(const T &f, const T &phase, T *x, size_t length)
{
    // generate signal
    for (size_t n = 0; n < length; ++n) {
        x[n] = dsp::sin(dsp::constants<float>::pi * dsp::load(f) * n +
                        dsp::load(phase));
    }
}

template <class Filter, typename T>
static void filterArray(const T &cut, T *x, size_t length)
{
    Filter filter;
    typename decltype(filter)::State filterstate;

    filter.setFreq(cut);
    filter.processBlock(dsp::Context(x, length), filterstate);
}

template <class Filter, typename T>
static void filterResArray(const T &cut, const T &res, T *x, size_t length)
{
    Filter filter;
    typename decltype(filter)::State filterstate;

    filter.setFreq(cut, res);
    filter.processBlock(dsp::Context(x, length), filterstate);
}

template <typename T> static auto rms(T *x, size_t length)
{
    auto rms = dsp::load(T(0));
    for (size_t n = 0; n < length; ++n) {
        rms += dsp::load(x[n]) * x[n];
    }
    rms = dsp::sqrt(rms / length);

    return rms;
}

TEST_CASE("VA Filters", "[dsp][va][onepole][lowpass]")
{

    // test filter reaction to DC signal
    SECTION("DC")
    {
        constexpr size_t Length = 1 << 10;
        using ft                = float;
        ft x[Length];

        auto cut = GENERATE(500.f / kNyquist, 1000.f / kNyquist,
                            2000.f / kNyquist, 15000.f / kNyquist);

#define TEST_DC(FILTER, TYPE, VAL)                                     \
    {                                                                  \
        int start  = 6.f / cut;                                        \
        int length = start + 20;                                       \
        generatesine(0.f, {dsp::constants<ft>::pi / 2.f}, x, length);  \
        filterArray<FILTER<ft, TYPE>>(cut, x, length);                 \
        for (auto n = start; n < length; ++n) {                        \
            REQUIRE_THAT(x[n], Catch::Matchers::WithinAbs(VAL, 1e-3)); \
        }                                                              \
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
        using ft                = float;
        constexpr size_t Length = 1 << 15;
        float cutGainm3db       = sqrtf(1.f / 2.f);
        float cutGainm12db      = 1 / 4.f;
        float power;
        ft x[Length];

        auto cut = GENERATE(take(6, random(0.f, 1.f)));

#define TEST_CUTOFF(FILTER, TYPE, GAIN)              \
    generatesine(cut, 0.f, x, Length);               \
    power = rms(x, Length);                          \
    filterArray<FILTER<ft, TYPE>>({cut}, x, Length); \
    REQUIRE_THAT(rms(x, Length) / power,             \
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
        constexpr size_t N      = 4;
        constexpr size_t Length = 1 << 10;
        using ft                = dsp::mfloat<N>;
        constexpr ft fs         = {0.97f, 0.964f, 0.9134f, 0.9273f};
        ft x[Length];
        ft powerOrig;
        ft power;

        auto cut = GENERATE(50.3f / kNyquist, 104.38f / kNyquist,
                            200.38f / kNyquist, 332.2f / kNyquist);

#define TEST_HIGHFREQ(FILTER, TYPE, VAL)                               \
    generatesine(fs, ft(0), x, Length);                                \
    powerOrig = rms(x, Length);                                        \
    filterArray<FILTER<ft, TYPE>>({cut, cut, cut, cut}, x, Length);    \
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
        using ft                = float;
        constexpr size_t Length = 1 << 15;
        ft power;
        ft x[Length];

        auto cut = GENERATE(take(4, random(0.05f, 0.5f)));
        auto R   = GENERATE(0.2f, take(4, random(0.08f, 2.f)));

        float cutGainm3db = 1.f / (2.f * R);

#define TEST_RES(FILTER, TYPE, GAIN)                         \
    generatesine(cut, 0.f, x, Length);                       \
    power = rms(x, Length);                                  \
    filterResArray<FILTER<ft, TYPE>>({cut}, {R}, x, Length); \
    REQUIRE_THAT(rms(x, Length) / power,                     \
                 Catch::Matchers::WithinAbs(GAIN, 5e-2));

        TEST_RES(dsp::va::SVF, dsp::va::kLowPass, cutGainm3db);
        TEST_RES(dsp::va::SVF, dsp::va::kHighPass, cutGainm3db);
        TEST_RES(dsp::va::SVF, dsp::va::kBandPass, 1.f);
    }
}
