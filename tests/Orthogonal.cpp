#include "dsp/Orthogonal.h"
#include "dsp/Context.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("Householder transform", "[dsp][householder]")
{
    using namespace Catch::Matchers;

    SECTION("N=4")
    {
        SECTION("float")
        {
            dsp::mfloat<4> x = {1.f, 2.f, 3.f, 4.f};
            auto y           = dsp::householder(dsp::load(x));

            float expected[] = {
                -8.f, // 1-2-3-4
                -6.f, //-1+2-3-4
                -4.f, //-1-2+3-4
                -2.f  //-1-2-3+4
            };

            CHECK(y[0] == expected[0] / 2.f);
            CHECK(y[1] == expected[1] / 2.f);
            CHECK(y[2] == expected[2] / 2.f);
            CHECK(y[3] == expected[3] / 2.f);
        }

#if DSP_MAX_VEC_SIZE / 8 > 2
        SECTION("double")
        {
            dsp::mdouble<4> x = {1., 2., 3., 4.};
            auto y            = dsp::householder(dsp::load(x));

            double expected[] = {
                -8., // 1-2-3-4
                -6., //-1+2-3-4
                -4., //-1-2+3-4
                -2.  //-1-2-3+4
            };

            CHECK(y[0] == expected[0] / 2.);
            CHECK(y[1] == expected[1] / 2.);
            CHECK(y[2] == expected[2] / 2.);
            CHECK(y[3] == expected[3] / 2.);
        }
#endif
    }
}

TEST_CASE("Hadamard transform", "[dsp][hadamard]")
{
    using namespace Catch::Matchers;

    SECTION("Identity case (L=1)")
    {
        float x = 5.f;
        auto y  = dsp::hadamard(dsp::load(x));
        CHECK(y == 5.f);
    }

    SECTION("Basic Hadamard case (N=2)")
    {
        SECTION("float")
        {
            dsp::mfloat<2> x      = {1.f, 2.f};
            auto y                = dsp::hadamard(dsp::load(x));
            constexpr auto isqrt2 = dsp::constants<float>::sqrt1_2;
            CHECK(y[0] == 3.f * isqrt2);  // 1 + 2
            CHECK(y[1] == -1.f * isqrt2); // 1 - 2
        }
#ifdef DSP_SIMD_DOUBLE
        SECTION("double")
        {
            dsp::mdouble<2> x     = {1., 2.};
            auto y                = dsp::hadamard(dsp::load(x));
            constexpr auto isqrt2 = dsp::constants<double>::sqrt1_2;
            CHECK(y[0] == 3. * isqrt2);  // 1 + 2
            CHECK(y[1] == -1. * isqrt2); // 1 - 2
        }
#endif
    }

    SECTION("N=4")
    {
        SECTION("float")
        {
            dsp::mfloat<4> x = {1.f, 2.f, 3.f, 4.f};
            auto y           = dsp::hadamard(dsp::load(x));

            float expected[] = {
                10.f, // 1+2+3+4
                -2.f, // 1-2+3-4
                -4.f, // 1+2-3-4
                0.f   // 1-2-3+4
            };

            CHECK(y[0] == expected[0] / 2.f);
            CHECK(y[1] == expected[1] / 2.f);
            CHECK(y[2] == expected[2] / 2.f);
            CHECK(y[3] == expected[3] / 2.f);
        }

#if DSP_MAX_VEC_SIZE / 4 > 4
        SECTION("double")
        {
            dsp::mdouble<4> x = {1., 2., 3., 4.};
            auto y            = dsp::hadamard(dsp::load(x));

            double expected[] = {
                10., // 1+2+3+4
                -2., // 1-2+3-4
                -4., // 1+2-3-4
                0.   // 1-2-3+4
            };

            CHECK(y[0] == expected[0] / 2.);
            CHECK(y[1] == expected[1] / 2.);
            CHECK(y[2] == expected[2] / 2.);
            CHECK(y[3] == expected[3] / 2.);
        }
#endif
    }

#if DSP_MAX_VEC_SIZE / 4 > 4
    SECTION("N=8")
    {
        dsp::mfloat<8> x = {1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f};
        const auto y     = hadamard(dsp::load(x));
        float norm       = sqrtf(1.f / 8.f);

        float expected[] = {
            36.f,  // 1+2+3+4+5+6+7+8
            -4.f,  // 1-2+3-4+5-6+7-8
            -8.f,  // 1+2-3-4+5+6-7-8
            0.f,   // 1-2-3+4+5-6-7+8
            -16.f, // 1+2+3+4-5-6-7-8
            0.f,   // 1-2+3-4-5+6-7+8
            0.f,   // 1+2-3-4-5-6+7+8
            0.f,   // 1-2-3+4-5+6+7-8
        };

        REQUIRE_THAT(y[0], WithinAbs(norm * expected[0], 1e-6));
        REQUIRE_THAT(y[1], WithinAbs(norm * expected[1], 1e-6));
        REQUIRE_THAT(y[2], WithinAbs(norm * expected[2], 1e-6));
        REQUIRE_THAT(y[3], WithinAbs(norm * expected[3], 1e-6));
        REQUIRE_THAT(y[4], WithinAbs(norm * expected[4], 1e-6));
        REQUIRE_THAT(y[5], WithinAbs(norm * expected[5], 1e-6));
        REQUIRE_THAT(y[6], WithinAbs(norm * expected[6], 1e-6));
        REQUIRE_THAT(y[7], WithinAbs(norm * expected[7], 1e-6));
    }
#endif
}

// TEST_CASE("Hadamard interpolation", "[dsp][hadamard][interpolation]")
//{
//     using namespace Catch::Matchers;
//     using namespace Catch::Generators;
//
//     SECTION("Identity t=0")
//     {
//         constexpr auto kN = 4;
//         auto h0           = dsp::hadamardInterpolMatrix<kN>(0.f);
//
//         REQUIRE_THAT(h0[0][0], WithinAbs(1.f, 1e-7));
//         REQUIRE_THAT(h0[0][1], WithinAbs(0.f, 1e-7));
//         REQUIRE_THAT(h0[0][2], WithinAbs(0.f, 1e-7));
//         REQUIRE_THAT(h0[0][3], WithinAbs(0.f, 1e-7));
//
//         REQUIRE_THAT(h0[1][0], WithinAbs(0.f, 1e-7));
//         REQUIRE_THAT(h0[1][1], WithinAbs(1.f, 1e-7));
//         REQUIRE_THAT(h0[1][2], WithinAbs(0.f, 1e-7));
//         REQUIRE_THAT(h0[1][3], WithinAbs(0.f, 1e-7));
//
//         REQUIRE_THAT(h0[2][0], WithinAbs(0.f, 1e-7));
//         REQUIRE_THAT(h0[2][1], WithinAbs(0.f, 1e-7));
//         REQUIRE_THAT(h0[2][2], WithinAbs(1.f, 1e-7));
//         REQUIRE_THAT(h0[2][3], WithinAbs(0.f, 1e-7));
//
//         REQUIRE_THAT(h0[3][0], WithinAbs(0.f, 1e-7));
//         REQUIRE_THAT(h0[3][1], WithinAbs(0.f, 1e-7));
//         REQUIRE_THAT(h0[3][2], WithinAbs(0.f, 1e-7));
//         REQUIRE_THAT(h0[3][3], WithinAbs(1.f, 1e-7));
//     }
//
//     SECTION("Hadamard matrix t=1")
//     {
//         constexpr auto kN = 4;
//         auto h1           = dsp::hadamardInterpolMatrix<kN>(1.f);
//
//         dsp::linalg::fMatrix<kN> h = {{{
//             {.5, .5, .5, .5},
//             {.5, -.5, .5, -.5},
//             {.5, .5, -.5, -.5},
//             {.5, -.5, -.5, .5},
//         }}};
//         REQUIRE_THAT(h1[0][0], WithinAbs(h[0][0], 1e-7));
//         REQUIRE_THAT(h1[0][1], WithinAbs(h[0][1], 1e-7));
//         REQUIRE_THAT(h1[0][2], WithinAbs(h[0][2], 1e-7));
//         REQUIRE_THAT(h1[0][3], WithinAbs(h[0][3], 1e-7));
//
//         REQUIRE_THAT(h1[1][0], WithinAbs(h[1][0], 1e-7));
//         REQUIRE_THAT(h1[1][1], WithinAbs(h[1][1], 1e-7));
//         REQUIRE_THAT(h1[1][2], WithinAbs(h[1][2], 1e-7));
//         REQUIRE_THAT(h1[1][3], WithinAbs(h[1][3], 1e-7));
//
//         REQUIRE_THAT(h1[2][0], WithinAbs(h[2][0], 1e-7));
//         REQUIRE_THAT(h1[2][1], WithinAbs(h[2][1], 1e-7));
//         REQUIRE_THAT(h1[2][2], WithinAbs(h[2][2], 1e-7));
//         REQUIRE_THAT(h1[2][3], WithinAbs(h[2][3], 1e-7));
//
//         REQUIRE_THAT(h1[3][0], WithinAbs(h[3][0], 1e-7));
//         REQUIRE_THAT(h1[3][1], WithinAbs(h[3][1], 1e-7));
//         REQUIRE_THAT(h1[3][2], WithinAbs(h[3][2], 1e-7));
//         REQUIRE_THAT(h1[3][3], WithinAbs(h[3][3], 1e-7));
//     }
//
//     SECTION("interpolation matrix is orthogonal")
//     {
//         constexpr auto kN = 4;
//         float t           = GENERATE(take(4, random(0.f, 1.f)));
//         auto h            = dsp::hadamardInterpolMatrix<kN>(t);
//
//         // H is ortogonal
//         REQUIRE_THAT(dot(h[0], h[0]), WithinAbs(1.f, 1e-6));
//         REQUIRE_THAT(dot(h[0], h[1]), WithinAbs(0.f, 1e-6));
//         REQUIRE_THAT(dot(h[0], h[2]), WithinAbs(0.f, 1e-6));
//         REQUIRE_THAT(dot(h[0], h[3]), WithinAbs(0.f, 1e-6));
//
//         REQUIRE_THAT(dot(h[1], h[1]), WithinAbs(1.f, 1e-6));
//         REQUIRE_THAT(dot(h[1], h[2]), WithinAbs(0.f, 1e-6));
//         REQUIRE_THAT(dot(h[1], h[3]), WithinAbs(0.f, 1e-6));
//
//         REQUIRE_THAT(dot(h[2], h[2]), WithinAbs(1.f, 1e-6));
//         REQUIRE_THAT(dot(h[2], h[3]), WithinAbs(0.f, 1e-6));
//
//         REQUIRE_THAT(dot(h[3], h[3]), WithinAbs(1.f, 1e-6));
//     }
// }
