#include "dsp/Noise.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>

using namespace Catch::Matchers;

template <typename T, int N = 1000> static void testNoise()
{
    using bt = dsp::baseType<T>;
    dsp::Noise<T> noise(3848230);

    auto sum = dsp::load(T{});
    auto var = sum;
    bt covar = 0;
    // run lfo
    for (int i = 0; i < N; ++i) {
        auto x = noise.process();
        sum += x;
        var += x * x;
        bt tmpcorr = 1;
        for (size_t i = 0; i < dsp::kTypeWidth<T>; ++i)
            tmpcorr *= dsp::get(x, i);
        covar += tmpcorr;
    }
    sum /= N;
    var /= N;
    covar /= N;

    for (size_t i = 0; i < dsp::kTypeWidth<T>; ++i) {
        REQUIRE_THAT(dsp::get(sum, i), WithinAbs(bt(0), bt(1e-1)));
        REQUIRE_THAT(dsp::get(var, i), WithinAbs(bt(1) / 3, bt(5e-2)));
    }
    REQUIRE_THAT(covar, WithinAbs(bt(0), bt(1e-1)));
}

TEST_CASE("Noise", "[dsp][noise]")
{
    testNoise<float>();
    testNoise<dsp::mfloat<2>>();
    testNoise<dsp::mdouble<2>>();
    testNoise<dsp::mfloat<4>>();
#if DSP_MAX_VEC_SIZE / 4 > 4
    testNoise<dsp::mfloat<8>>();
#endif
}
