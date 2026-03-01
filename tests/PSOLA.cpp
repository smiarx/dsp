#include <catch2/catch_all.hpp>
#include <dsp/delay/Delay.hpp>
#include <dsp/pitch/PSOLA.hpp>
#include <iostream>

using namespace Catch::Matchers;

TEST_CASE("PSOLA")
{
    // dsp::PSOLA<float, 8> psola{};
    // dsp::DelayLine<200> delaylineAnalysis;
    // dsp::DelayLine<500, nextTo(delaylineAnalysis)> delaylineSynth;
    // dsp::Buffer<float, nextTo(delaylineSynth)> buffer;
    // std::array<float, decltype(buffer)::kSize> bufferData{};

    // buffer.setData(bufferData.data());

    // constexpr auto kN            = 1000;
    // constexpr auto kPeriod       = 80.f;
    // constexpr auto kHarmonicGain = 0.7f;
    // constexpr auto kShift        = 0.3f;

    ////auto startId = buffer.getId();

    // for (int i = 0; i < kN; ++i) {
    //     auto x = dsp::cos(2 * dsp::constants<float>::pi / kPeriod *
    //                       static_cast<float>(i)) +
    //              kHarmonicGain * dsp::cos(2 * dsp::constants<float>::pi * 2.f
    //              /
    //                                       kPeriod * static_cast<float>(i));
    //     dsp::BufferContext ctxt(&x, 1, buffer);

    //    if (i == 85) {
    //        psola.start(ctxt, kPeriod);

    //        //auto marker = psola.getMarker(0);
    //        //REQUIRE_THAT(
    //        //    marker.position_,
    //        //    WithinAbs(static_cast<float>(startId) + kPeriod, 1e-10));
    //        //REQUIRE_THAT(marker.period_, WithinAbs(kPeriod, 1e-10));

    //        //auto fade = psola.getFade();
    //        REQUIRE((psola.getState() == psola.kFadeIn));
    //        //REQUIRE_THAT(fade.rate_, WithinAbs(1.f / kPeriod, 1e-10));

    //        //auto grain = psola.getGrain(0);
    //        REQUIRE((psola.getNumGrains() == 1));
    //        //REQUIRE_THAT(grain.rate_, WithinAbs(1.f / kPeriod, 1e-10));
    //        //REQUIRE_THAT(
    //        //    grain.phase_,
    //        //    WithinAbs(-1.f + (ctxt.idDiff(marker.position_)) / kPeriod,
    //        //              1e-10));
    //    } else if (i == 820) {
    //        psola.stop(ctxt);
    //    }

    //    psola.setFreq(2.f/kPeriod * kShift);
    //    psola.process(ctxt, delaylineAnalysis, delaylineSynth);
    //    buffer.nextBlock(ctxt);

    //    //std::cout << ctxt.getInput() << "\n";
    //}
}
