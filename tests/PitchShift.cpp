#include "dsp/PitchShift.h"
#include "dsp/Buffer.h"

#include <catch2/catch_all.hpp>
#include <iostream>

TEST_CASE("Pitch Shift", "[dsp][pitchshift]")
{

    using namespace Catch::Matchers;

    /*
     * pitch shift sines and check if result is correct pitch
     */
    SECTION("shift sine")
    {

        dsp::PitchShift<float> ps;
        dsp::DelayLine<3000> dl;
        dsp::Buffer<float, nextTo(dl)> buffer;
        float bufferData[decltype(buffer)::kSize]{};
        buffer.setData(bufferData);

        constexpr auto kN = 100;
        float x[kN]{};

        auto freq = 0.05f;
        ps.setFreq(freq);

        dsp::BufferContext ctxt(x, kN, buffer);

        for (size_t n = 0; n < kN; ++n)
            x[n] = dsp::sin(dsp::constants<float>::pi * freq *
                            static_cast<float>(n));

        auto shift = GENERATE(0.364f, 0.564f, 0.738f, 1.f, 1.232, 1.564f, 2.f);
        ps.setShift(shift);
        int n = 0;
        CTXTRUN(ctxt)
        {
            ps.process(ctxt, dl);
            // expected value
            auto e = dsp::sin(dsp::constants<float>::pi * freq * shift *
                              (static_cast<float>(n) - 3.f / shift + 1));

            // first few samples wont be right
            if (n > 50) REQUIRE_THAT(ctxt.getInput(), WithinAbs(e, 1e-3));
            ++n;
        };
    }
}
