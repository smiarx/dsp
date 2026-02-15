#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>

#include "dsp/Buffer.h"
#include "dsp/MultiRate.h"

TEST_CASE("MultiRate", "[dsp][multirate]")
{
    constexpr auto kN     = 400;
    constexpr auto kOrder = 24;

    using ft = float;
    dsp::MultiRate<ft, kOrder, 6> multirate;

    decltype(multirate)::DLDecimate<0> decState;
    decltype(multirate)::DLInterpolate<0> interpState;

    // delayline buffers
    dsp::Buffer<ft, nextTo(interpState)> bufDec{};
    std::array<ft, decltype(bufDec)::kSize> bufDecData{};
    bufDec.setData(bufDecData.data());

    dsp::Buffer<ft, nextTo(decState)> buffer{};
    std::array<ft, decltype(buffer)::kSize> bufferData{};
    buffer.setData(bufferData.data());

    std::array<ft, kN> x{};
    std::array<ft, kN> xDec{};
    std::array<ft, kN> xInterp{};

    // low frequency signal
    for (size_t i = 0; i < kN; ++i) {
        x[i] = dsp::sin(2 * dsp::constants<ft>::pi * ft(200. / 48000) * ft(i));
    }

    auto ctxt       = dsp::BufferContext(x.data(), kN, buffer);
    auto ctxtDec    = dsp::BufferContext(xDec.data(), kN, bufDec);
    auto ctxtInterp = dsp::Context(xInterp.data(), kN);

    auto kRate = GENERATE(2u, 3u, 4u, 5u, 6u);
    multirate.decimate(kRate, ctxt, ctxtDec, decState, 0);
    multirate.interpolate(kRate, ctxtDec, ctxtInterp, interpState, 0);

    // decimate & interpolate will reproduce same signal delayed
    auto kDelay = static_cast<size_t>(multirate.getDelay(kRate));
    for (size_t i = kDelay + 100; i < kN; ++i) {
        REQUIRE_THAT(xInterp[i],
                     Catch::Matchers::WithinAbs(x[i - kDelay], ft(1e-6)));
    }
}
