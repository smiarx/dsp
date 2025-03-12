#include "dsp/FastMath.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>


TEST_CASE("Walsh-Hadamard", "[dsp][math][hadamard]")
{
using namespace Catch::Matchers;

    SECTION("Identity case (L=1)")
    {
        dsp::fSample<1> x = {5.f};
        auto y            = dsp::hadamard(x);
        CHECK(y[0] == 5.f);
    }

    SECTION("Basic Hadamard case (N=2)")
    {
        dsp::fSample<2> x = {1.f, 2.f};
        auto y            = hadamard(x);
        float isqrt2 = std::sqrt(1.f/2.f);
        CHECK(y[0] == 3.f*isqrt2);  // 1 + 2
        CHECK(y[1] == -1.f*isqrt2); // 1 - 2
    }

    SECTION("N=4")
    {
        dsp::fSample<4> x = {1.f, 2.f, 3.f, 4.f};
        auto y            = hadamard(x);

        dsp::fSample<4> expected = {
            10.f, // 1+2+3+4
            -2.f, // 1-2+3-4
            -4.f, // 1+2-3-4
            0.f   // 1-2-3+4
        };

        CHECK(y[0] == expected[0]/2.f);
        CHECK(y[1] == expected[1]/2.f);
        CHECK(y[2] == expected[2]/2.f);
        CHECK(y[3] == expected[3]/2.f);
    }

#if defined(__AVX)
    SECTION("N=8")
    {
        dsp::fSample<8> x = {1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f};
        auto y            = hadamard(x);
        float norm = sqrtf(1.f/8.f);

        dsp::fSample<8> expected = {
            36.f,  // 1+2+3+4+5+6+7+8
            -4.f,  // 1-2+3-4+5-6+7-8
            -8.f,  // 1+2-3-4+5+6-7-8
            0.f,   // 1-2-3+4+5-6-7+8
            -16.f, // 1+2+3+4-5-6-7-8
            0.f,   // 1-2+3-4-5+6-7+8
            0.f,   // 1+2-3-4-5-6+7+8
            0.f,   // 1-2-3+4-5+6+7-8
        };

        REQUIRE_THAT(y[0] ,WithinAbs(norm*expected[0], 1e-6));
        REQUIRE_THAT(y[1] ,WithinAbs(norm*expected[1], 1e-6));
        REQUIRE_THAT(y[2] ,WithinAbs(norm*expected[2], 1e-6));
        REQUIRE_THAT(y[3] ,WithinAbs(norm*expected[3], 1e-6));
        REQUIRE_THAT(y[4] ,WithinAbs(norm*expected[4], 1e-6));
        REQUIRE_THAT(y[5] ,WithinAbs(norm*expected[5], 1e-6));
        REQUIRE_THAT(y[6] ,WithinAbs(norm*expected[6], 1e-6));
        REQUIRE_THAT(y[7] ,WithinAbs(norm*expected[7], 1e-6));
    }
#endif
}
