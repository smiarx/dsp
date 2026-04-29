#include <catch2/catch_all.hpp>
#include <dsp/music/Pitch.hpp>

using namespace dsp;
using namespace dsp::units;

TEST_CASE("Pitch", "[music][pitch]")
{
    using namespace Catch::Matchers;

    auto ref = PitchReference{440._Hz};

    auto pitch = ref.toPitch(440._Hz);
    REQUIRE_THAT(pitch.get(), WithinAbs(0., 1e-9));

    auto pitch2 = ref.toPitch(880._Hz);
    REQUIRE_THAT(pitch2.get(),
                 WithinAbs(Pitch<double>::octaveUnit().get(), 1e-9));

    auto pitch_2 = ref.toPitch(220._Hz);
    REQUIRE_THAT(pitch_2.get(),
                 WithinAbs(-Pitch<double>::octaveUnit().get(), 1e-9));

    auto freqSqrt2 = ref.toFreq(Pitch<double>::octaveUnit() * 0.5);
    REQUIRE_THAT(freqSqrt2.get(), WithinAbs(440. * std::sqrt(2.), 1e-9));
}
