#include "dsp/LFO.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace dsp::loadfuncs;
using namespace Catch::Matchers;

template <typename T, int N = 200>
static void testSine(const T &freq, const T &phase = {})
{
    dsp::lfo::Sine<T> lfo(freq, phase);

    auto f = load(freq);
    auto p = load(phase);
    // run lfo
    for (int i = 0; i < N - 100; ++i) {
        auto x = lfo.process();
        (void)x;
    }
    for (int i = N - 100; i < N; ++i) {
        auto x = lfo.process();
        auto expect =
            dsp::sin(dsp::constants<T>::pi * (f * static_cast<T>(i) + p * 2));
        for (size_t i = 0; i < dsp::kTypeWidth<T>; ++i)
            REQUIRE_THAT(dsp::get(x, i), WithinAbs(dsp::get(expect, i), 1e-3f));
    }
}

TEST_CASE("LFO Sine", "[dsp][lfo][sine]")
{
    SECTION("float") { testSine(0.001f); }
    SECTION("float with phase") { testSine(0.005f, 0.684f); }
    SECTION("float x 2") { testSine<dsp::mfloat<2>>({0.002, 0.003}); }
    SECTION("float x 2 with phase")
    {
        testSine<dsp::mfloat<2>>({0.0002339, 0.000332}, {0.3982, 0.9883});
    }
    SECTION("double x 2") { testSine<dsp::mdouble<2>>({0.000212, 0.0008}); }
    SECTION("float x 4")
    {
        testSine<dsp::mfloat<4>>({0.00143, 0.0008, 0.0032, 0.0011});
    }

    SECTION("long run") { testSine<double, 48000>(0.000354); }
}

template <typename T, int N = 200> static void testParabolic(const T &freq)
{
    dsp::lfo::Parabolic<T> lfo(freq);

    auto f = load(freq);
    // run lfo
    for (int i = 0; i < 100; ++i) {
        auto x = lfo.process();

        T expect     = i * f + 0.5;
        auto iexpect = dsp::toInt(expect);
        expect -= iexpect + 0.5;

        expect = 4 * (expect);
        expect = expect * (2 - dsp::abs(expect));

        for (size_t i = 0; i < dsp::kTypeWidth<T>; ++i)
            REQUIRE_THAT(dsp::get(x, i), WithinAbs(dsp::get(expect, i), 1e-3f));
    }
}

TEST_CASE("LFO Parabolic", "[dsp][lfo][parabolic]")
{
    testParabolic<float>(0.0004f);
    testParabolic<double>(0.005f);
}
