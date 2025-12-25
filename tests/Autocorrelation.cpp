#include "dsp/Autocorrelation.h"
#include "dsp/Buffer.h"
#include "dsp/FastMath.h"
#include <catch2/catch_all.hpp>

TEST_CASE("Autocorrelation")
{
    dsp::NSDF<float, 100> nsdf;
    decltype(nsdf)::DL<0> dl;
    dsp::Buffer<float, nextTo(dl)> buffer;
    float bufferData[decltype(buffer)::kSize]{};
    buffer.setData(bufferData);

    constexpr auto kN = 300;
    float x[kN]{};
    int period = GENERATE(5, 18, 30, 60);

    for (size_t n = 0; n < kN; ++n) {
        x[n] = dsp::sin(2 * dsp::constants<float>::pi /
                        static_cast<float>(period) * static_cast<float>(n));
    }

    dsp::BufferContext ctxt(x, kN, buffer);
    CTXTRUN(ctxt) { nsdf.process(ctxt, dl); };

    // multiples of period have nsdf of 1
    for (size_t t = period; t < decltype(nsdf)::kTauLength; t += period) {
        REQUIRE_THAT(nsdf.get(t), Catch::Matchers::WithinAbs(1.f, 1e-6f));
    }
}
