#include "dsp/LinAlg.h"
#include "dsp/Orthogonal.h"
#include "dsp/Signal.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstddef>

template <size_t N> static float dot(dsp::fData<N> x, dsp::fData<N> y)
{
    float val{};
    for (size_t i = 0; i < N; ++i) {
        val += x[i] * y[i];
    }
    return val;
}

TEST_CASE("Hadamard transform", "[dsp][hadamard]")
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
        float isqrt2      = std::sqrt(1.f / 2.f);
        CHECK(y[0] == 3.f * isqrt2);  // 1 + 2
        CHECK(y[1] == -1.f * isqrt2); // 1 - 2
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

        CHECK(y[0] == expected[0] / 2.f);
        CHECK(y[1] == expected[1] / 2.f);
        CHECK(y[2] == expected[2] / 2.f);
        CHECK(y[3] == expected[3] / 2.f);
    }

#if defined(__AVX)
    SECTION("N=8")
    {
        dsp::fSample<8> x = {1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f};
        auto y            = hadamard(x);
        float norm        = sqrtf(1.f / 8.f);

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

TEST_CASE("Hadamard interpolation", "[dsp][hadamard][interpolation]")
{
    using namespace Catch::Matchers;
    using namespace Catch::Generators;

    SECTION("Identity t=0")
    {
        constexpr auto kN = 4;
        auto h0           = dsp::hadamardInterpolMatrix<kN>(0.f);

        REQUIRE_THAT(h0[0][0], WithinAbs(1.f, 1e-7));
        REQUIRE_THAT(h0[0][1], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(h0[0][2], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(h0[0][3], WithinAbs(0.f, 1e-7));

        REQUIRE_THAT(h0[1][0], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(h0[1][1], WithinAbs(1.f, 1e-7));
        REQUIRE_THAT(h0[1][2], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(h0[1][3], WithinAbs(0.f, 1e-7));

        REQUIRE_THAT(h0[2][0], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(h0[2][1], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(h0[2][2], WithinAbs(1.f, 1e-7));
        REQUIRE_THAT(h0[2][3], WithinAbs(0.f, 1e-7));

        REQUIRE_THAT(h0[3][0], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(h0[3][1], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(h0[3][2], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(h0[3][3], WithinAbs(1.f, 1e-7));
    }

    SECTION("Hadamard matrix t=1")
    {
        constexpr auto kN = 4;
        auto h1           = dsp::hadamardInterpolMatrix<kN>(1.f);

        dsp::linalg::fMatrix<kN> h = {{{
            {.5, .5, .5, .5},
            {.5, -.5, .5, -.5},
            {.5, .5, -.5, -.5},
            {.5, -.5, -.5, .5},
        }}};
        REQUIRE_THAT(h1[0][0], WithinAbs(h[0][0], 1e-7));
        REQUIRE_THAT(h1[0][1], WithinAbs(h[0][1], 1e-7));
        REQUIRE_THAT(h1[0][2], WithinAbs(h[0][2], 1e-7));
        REQUIRE_THAT(h1[0][3], WithinAbs(h[0][3], 1e-7));

        REQUIRE_THAT(h1[1][0], WithinAbs(h[1][0], 1e-7));
        REQUIRE_THAT(h1[1][1], WithinAbs(h[1][1], 1e-7));
        REQUIRE_THAT(h1[1][2], WithinAbs(h[1][2], 1e-7));
        REQUIRE_THAT(h1[1][3], WithinAbs(h[1][3], 1e-7));

        REQUIRE_THAT(h1[2][0], WithinAbs(h[2][0], 1e-7));
        REQUIRE_THAT(h1[2][1], WithinAbs(h[2][1], 1e-7));
        REQUIRE_THAT(h1[2][2], WithinAbs(h[2][2], 1e-7));
        REQUIRE_THAT(h1[2][3], WithinAbs(h[2][3], 1e-7));

        REQUIRE_THAT(h1[3][0], WithinAbs(h[3][0], 1e-7));
        REQUIRE_THAT(h1[3][1], WithinAbs(h[3][1], 1e-7));
        REQUIRE_THAT(h1[3][2], WithinAbs(h[3][2], 1e-7));
        REQUIRE_THAT(h1[3][3], WithinAbs(h[3][3], 1e-7));
    }

    SECTION("interpolation matrix is orthogonal")
    {
        constexpr auto kN = 4;
        float t           = GENERATE(take(4, random(0.f, 1.f)));
        auto h            = dsp::hadamardInterpolMatrix<kN>(t);

        // H is ortogonal
        REQUIRE_THAT(dot(h[0], h[0]), WithinAbs(1.f, 1e-6));
        REQUIRE_THAT(dot(h[0], h[1]), WithinAbs(0.f, 1e-6));
        REQUIRE_THAT(dot(h[0], h[2]), WithinAbs(0.f, 1e-6));
        REQUIRE_THAT(dot(h[0], h[3]), WithinAbs(0.f, 1e-6));

        REQUIRE_THAT(dot(h[1], h[1]), WithinAbs(1.f, 1e-6));
        REQUIRE_THAT(dot(h[1], h[2]), WithinAbs(0.f, 1e-6));
        REQUIRE_THAT(dot(h[1], h[3]), WithinAbs(0.f, 1e-6));

        REQUIRE_THAT(dot(h[2], h[2]), WithinAbs(1.f, 1e-6));
        REQUIRE_THAT(dot(h[2], h[3]), WithinAbs(0.f, 1e-6));

        REQUIRE_THAT(dot(h[3], h[3]), WithinAbs(1.f, 1e-6));
    }
}
