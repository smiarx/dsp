#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>
#include <dsp/core/Context.hpp>
#include <dsp/core/SampleRate.hpp>
#include <dsp/filters/VAFilters.hpp>
#include <dsp/math/FastMath.hpp>

using namespace dsp::units;

static constexpr dsp::SampleRate<float> kSR{48000._Hzf};

template <typename T, dsp::freq::SampleRateIndependent Freq>
static void generatesine(const Freq &f, const T &phase, T *x, size_t length)
{
    auto angular = dsp::freq::Angular{f};
    // generate signal
    for (size_t n = 0; n < length; ++n) {
        x[n] = dsp::sin(dsp::load(static_cast<T>(angular.get())) * T(n) +
                        dsp::load(phase));
    }
}

template <class Filter, typename T, dsp::freq::SampleRateIndependent Freq>
static void filterArray(const Freq &cut, T *x, size_t length)
{
    Filter filter;
    typename decltype(filter)::State filterstate;

    filter.setFreq(cut);
    filter.processBlock(dsp::Context(x, int(length)), filterstate);
}

template <class Filter, typename T, typename U>
static void filterResArray(const dsp::freq::Angular<U> &cut, const T &res, T *x,
                           size_t length)
{
    Filter filter;
    typename decltype(filter)::State filterstate;

    filter.setFreq(cut, res);
    filter.processBlock(dsp::Context(x, int(length)), filterstate);
}

template <typename T> static auto rms(T *x, size_t length)
{
    auto rms = dsp::load(T(0));
    for (size_t n = 0; n < length; ++n) {
        rms += dsp::load(x[n]) * x[n];
    }
    rms = dsp::sqrt(rms / T(length));

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

        auto cut = GENERATE(500._Hzf, 1000._Hzf, 2000._Hzf, 15000._Hzf);

        auto angular = kSR.hzToAngular(cut);
        auto start   = size_t(30.f / angular.get());
        auto length  = start + 20;
        generatesine(0._nyquistf, {dsp::constants<ft>::pi / 2.f}, x, length);

#define TEST_DC(FILTER, TYPE, VAL)                                      \
    {                                                                   \
        filterArray<FILTER<ft, TYPE>>(kSR.hzToAngular(cut), x, length); \
        for (size_t n = start; n < length; ++n) {                       \
            REQUIRE_THAT(x[n], Catch::Matchers::WithinAbs(VAL, 1e-3));  \
        }                                                               \
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

        auto cut =
            dsp::freq::NyquistRatio(GENERATE(take(6, random(0.f, 0.95f))));

#define TEST_CUTOFF(FILTER, TYPE, GAIN)            \
    generatesine(cut, 0.f, x, Length);             \
    power = rms(x, Length);                        \
    filterArray<FILTER<ft, TYPE>>(cut, x, Length); \
    REQUIRE_THAT(rms(x, Length) / power,           \
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
        auto fs =
            dsp::freq::NyquistRatio<ft>{{0.97f, 0.964f, 0.9134f, 0.9273f}};
        ft x[Length];
        ft powerOrig;
        ft power;

        auto cut = GENERATE(50.3_Hz, 104.38_Hz, 200.38_Hz, 332.2_Hz);

#define TEST_HIGHFREQ(FILTER, TYPE, VAL)                               \
    generatesine(fs, ft(0), x, Length);                                \
    powerOrig = rms(x, Length);                                        \
    filterArray<FILTER<ft, TYPE>>(kSR.hzToAngular(cut), x, Length);    \
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

        auto cut =
            dsp::freq::Angular<ft>{GENERATE(take(4, random(0.05f, 0.5f)))};
        auto R = GENERATE(0.2f, take(4, random(0.08f, 2.f)));

        float cutGainm3db = 1.f / (2.f * R);

#define TEST_RES(FILTER, TYPE, GAIN)                       \
    generatesine(cut, 0.f, x, Length);                     \
    power = rms(x, Length);                                \
    filterResArray<FILTER<ft, TYPE>>(cut, {R}, x, Length); \
    REQUIRE_THAT(rms(x, Length) / power,                   \
                 Catch::Matchers::WithinAbs(GAIN, 5e-2));

        TEST_RES(dsp::va::SVF, dsp::va::kLowPass, cutGainm3db);
        TEST_RES(dsp::va::SVF, dsp::va::kHighPass, cutGainm3db);
        TEST_RES(dsp::va::SVF, dsp::va::kBandPass, 1.f);
    }
}
