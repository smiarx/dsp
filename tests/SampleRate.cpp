#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <dsp/core/SampleRate.hpp>

using Catch::Matchers::WithinRel;

using Hz      = dsp::freq::Hz<double>;
using Rad     = dsp::freq::Angular<double>;
using Nyquist = dsp::freq::NyquistRatio<double>;
using Seconds = dsp::time::Seconds<double>;
using Samples = dsp::time::Samples<double>;

TEST_CASE("SampleRate construction and get()", "[samplerate][basic]")
{
    // Arrange
    Hz srHz{48000.0};

    // Act
    dsp::SampleRate<double> sr{srHz};

    // Assert
    REQUIRE(sr.get().get() == 48000.0);
}

TEST_CASE("SampleRate set()", "[samplerate][basic]")
{
    dsp::SampleRate<double> sr{Hz{44100.0}};

    sr.set(Hz{48000.0});

    REQUIRE(sr.get().get() == 48000.0);
}

TEST_CASE("Hz to Angular conversion", "[samplerate][conversion]")
{
    dsp::SampleRate<double> sr{Hz{48000.0}};
    Hz freq{12000.0}; // quarter of sample rate

    auto ang = sr.hzToAngular(freq);

    // expected: pi * (2/sr) * freq = pi * (freq / (sr/2))
    double expected = dsp::constants<double>::pi * (2.0 / 48000.0) * 12000.0;

    REQUIRE_THAT(ang.get(), WithinRel(expected));
}

TEST_CASE("Hz to NyquistRatio conversion", "[samplerate][conversion]")
{
    dsp::SampleRate<double> sr{Hz{48000.0}};
    Hz freq{12000.0};

    auto ratio = sr.hzToNyquistRatio(freq);

    // expected: freq / (sr/2)
    double expected = 2.0 * 12000.0 / 48000.0;

    REQUIRE_THAT(ratio.get(), WithinRel(expected));
}

TEST_CASE("Seconds to Samples conversion", "[samplerate][conversion]")
{
    dsp::SampleRate<double> sr{Hz{48000.0}};
    Seconds sec{0.5};

    auto samples = sr.secsToSamples(sec);

    REQUIRE_THAT(samples.get(), WithinRel(24000.0));
}

TEST_CASE("Hz to period (samples)", "[samplerate][conversion]")
{
    dsp::SampleRate<double> sr{Hz{48000.0}};
    Hz freq{1000.0};

    auto period = sr.hzToPeriod(freq);

    REQUIRE_THAT(period.get(), WithinRel(48.0));
}

TEST_CASE("Period to Hz conversion (API sanity check)",
          "[samplerate][conversion]")
{
    dsp::SampleRate<double> sr{Hz{48000.0}};
    Samples period{48.0};

    auto result = sr.periodToHz(period);

    STATIC_REQUIRE(std::is_same_v<decltype(result), Hz>);

    // Numerically it's Hz, but wrapped incorrectly
    REQUIRE_THAT(result.get(), WithinRel(1000.0));
}

TEST_CASE("SampleRate scaling operators", "[samplerate][operators]")
{
    dsp::SampleRate<double> sr{Hz{48000.0}};

    SECTION("multiply")
    {
        auto scaled = sr * 2.0;

        REQUIRE_THAT(scaled.get().get(), WithinRel(96000.0));
    }

    SECTION("divide")
    {
        auto scaled = sr / 2.0;

        REQUIRE_THAT(scaled.get().get(), WithinRel(24000.0));
    }
}

TEST_CASE("SampleRate ratio operator", "[samplerate][operators]")
{
    dsp::SampleRate<double> a{Hz{48000.0}};
    dsp::SampleRate<double> b{Hz{24000.0}};

    auto ratio = a / b;

    REQUIRE_THAT(ratio, WithinRel(2.0));
}

TEST_CASE("constexpr behavior", "[samplerate][constexpr]")
{
    constexpr dsp::SampleRate<double> sr{Hz{48000.0}};
    constexpr auto ang = sr.hzToAngular(Hz{12000.0});

    STATIC_REQUIRE(ang.get() > 0.0);
}
