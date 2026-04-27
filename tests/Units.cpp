#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <dsp/core/Units.hpp>

using Catch::Approx;
using Catch::Matchers::WithinRel;

using Hz      = dsp::freq::Hz<double>;
using Rad     = dsp::freq::Angular<double>;
using Nyquist = dsp::freq::NyquistRatio<double>;
using Samples = dsp::time::Samples<double>;

TEST_CASE("Unit stores and returns value", "[unit][basic]")
{
    Hz hz{440.0};
    auto value = hz.get();
    REQUIRE(value == 440.0);
}

TEST_CASE("Scalar operations preserve unit type", "[unit][arithmetic]")
{
    Hz hz{100.0};

    SECTION("multiplication")
    {
        auto result = hz * 2.0;

        REQUIRE(result.get() == 200.0);
        STATIC_REQUIRE(std::is_same_v<decltype(result), Hz>);
    }

    SECTION("division")
    {
        auto result = hz / 2.0;

        REQUIRE(result.get() == 50.0);
        STATIC_REQUIRE(std::is_same_v<decltype(result), Hz>);
    }

    SECTION("compound multiplication")
    {
        hz *= 3.0;
        REQUIRE(hz.get() == 300.0);
    }
}

TEST_CASE("Scalar * Unit works symmetrically", "[unit][operators]")
{
    Hz hz{10.0};

    auto result = 3.0 * hz;

    REQUIRE(result.get() == 30.0);
    STATIC_REQUIRE(std::is_same_v<decltype(result), Hz>);
}

TEST_CASE("Unit addition and subtraction", "[unit][arithmetic]")
{
    Hz a{100.0};
    Hz b{40.0};

    SECTION("addition")
    {
        auto result = a + b;
        REQUIRE(result.get() == 140.0);
    }

    SECTION("subtraction")
    {
        auto result = a - b;
        REQUIRE(result.get() == 60.0);
    }
}

TEST_CASE("Unit divided by unit returns scalar", "[unit][arithmetic]")
{
    Hz a{100.0};
    Hz b{20.0};

    auto result = a / b;

    REQUIRE(result == 5.0);
    STATIC_REQUIRE(std::is_same_v<decltype(result), double>);
}

TEST_CASE("Equality operators", "[unit][comparison]")
{
    Hz a{100.0};
    Hz b{100.0};
    Hz c{50.0};

    REQUIRE(a == b);
    REQUIRE(a != c);
}

TEST_CASE("Cross-type operations deduce correct type", "[unit][types]")
{
    dsp::freq::Hz<float> a{100.0f};
    dsp::freq::Hz<double> b{50.0};

    auto result = a + b;

    REQUIRE(result.get() == Approx(150.0));
    STATIC_REQUIRE(std::is_same_v<decltype(result.get()), double>);
}

TEST_CASE("min/max on units", "[unit][algorithms]")
{
    Hz a{10.0};
    Hz b{20.0};

    SECTION("min")
    {
        auto result = dsp::min(a, b);
        REQUIRE(result.get() == 10.0);
    }

    SECTION("max")
    {
        auto result = dsp::max(a, b);
        REQUIRE(result.get() == 20.0);
    }
}

TEST_CASE("Angular <-> NyquistRatio conversions", "[unit][conversion]")
{
    Nyquist ratio{0.5};

    SECTION("Nyquist -> Angular")
    {
        Rad ang{ratio};

        REQUIRE_THAT(ang.get(), WithinRel(dsp::constants<double>::pi * 0.5));
    }

    SECTION("Angular -> Nyquist")
    {
        Rad ang{dsp::constants<double>::pi};

        Nyquist r{ang};

        REQUIRE_THAT(r.get(), WithinRel(1.0));
    }
}

TEST_CASE("NyquistRatio to Samples conversion", "[unit][conversion]")
{
    Nyquist r{0.5};

    auto samples = r.toPeriod();

    REQUIRE_THAT(samples.get(), WithinRel(4.0)); // 2 / 0.5
}

TEST_CASE("Samples to NyquistRatio conversion", "[unit][conversion]")
{
    Samples s{4.0};

    auto ratio = s.toNyquistRatio();

    REQUIRE_THAT(ratio.get(), WithinRel(0.5));
}

TEST_CASE("User-defined literals", "[unit][literals]")
{
    using namespace dsp::units;

    auto hz  = 440.0_Hz;
    auto rad = 3.14159_rad;
    auto ny  = 0.5_nyquist;
    auto sec = 1.0_s;
    auto ms  = 1.0_ms;

    REQUIRE(hz.get() == 440.0);
    REQUIRE_THAT(rad.get(), WithinRel(3.14159));
    REQUIRE(ny.get() == 0.5);
    REQUIRE(sec.get() == 1.0);
    REQUIRE_THAT(ms.get(), WithinRel(0.001));
}

TEST_CASE("Warp functions", "[unit][warp]")
{
    Nyquist r{0.5};

    auto result = dsp::freq::warp(r);

    REQUIRE_THAT(result,
                 WithinRel(std::tan(dsp::constants<double>::pi / 2 * 0.5)));
}

TEST_CASE("constexpr behavior", "[unit][constexpr]")
{
    constexpr Hz kHz{100.0};
    constexpr auto kDoubled = kHz * 2.0;

    STATIC_REQUIRE(kDoubled.get() == 200.0);
}
