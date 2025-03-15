#include "dsp/Hadamard.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

template <size_t N> float dot(dsp::fData<N> x, dsp::fData<N> y)
{
    float val{};
    for (size_t i = 0; i < N; ++i) {
        val += x[i] * y[i];
    }
    return val;
}

TEST_CASE("Hadamard", "[dsp][hadamard]")
{
    using namespace Catch::Matchers;
    using namespace Catch::Generators;

    SECTION("interpolation matrix t=0")
    {
        constexpr auto N = 4;
        auto H0          = dsp::hadamardInterpolMatrix<N>(0.f);

        REQUIRE_THAT(H0[0][0], WithinAbs(1.f, 1e-7));
        REQUIRE_THAT(H0[0][1], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(H0[0][2], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(H0[0][3], WithinAbs(0.f, 1e-7));

        REQUIRE_THAT(H0[1][0], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(H0[1][1], WithinAbs(1.f, 1e-7));
        REQUIRE_THAT(H0[1][2], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(H0[1][3], WithinAbs(0.f, 1e-7));

        REQUIRE_THAT(H0[2][0], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(H0[2][1], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(H0[2][2], WithinAbs(1.f, 1e-7));
        REQUIRE_THAT(H0[2][3], WithinAbs(0.f, 1e-7));

        REQUIRE_THAT(H0[3][0], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(H0[3][1], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(H0[3][2], WithinAbs(0.f, 1e-7));
        REQUIRE_THAT(H0[3][3], WithinAbs(1.f, 1e-7));
    }

    SECTION("interpolation matrix t=1")
    {
        constexpr auto N = 4;
        auto H1          = dsp::hadamardInterpolMatrix<N>(1.f);

        dsp::linalg::fMatrix<N> H = {{{
            {.5, .5, .5, .5},
            {.5, -.5, .5, -.5},
            {.5, .5, -.5, -.5},
            {.5, -.5, -.5, .5},
        }}};
        REQUIRE_THAT(H1[0][0], WithinAbs(H[0][0], 1e-7));
        REQUIRE_THAT(H1[0][1], WithinAbs(H[0][1], 1e-7));
        REQUIRE_THAT(H1[0][2], WithinAbs(H[0][2], 1e-7));
        REQUIRE_THAT(H1[0][3], WithinAbs(H[0][3], 1e-7));

        REQUIRE_THAT(H1[1][0], WithinAbs(H[1][0], 1e-7));
        REQUIRE_THAT(H1[1][1], WithinAbs(H[1][1], 1e-7));
        REQUIRE_THAT(H1[1][2], WithinAbs(H[1][2], 1e-7));
        REQUIRE_THAT(H1[1][3], WithinAbs(H[1][3], 1e-7));

        REQUIRE_THAT(H1[2][0], WithinAbs(H[2][0], 1e-7));
        REQUIRE_THAT(H1[2][1], WithinAbs(H[2][1], 1e-7));
        REQUIRE_THAT(H1[2][2], WithinAbs(H[2][2], 1e-7));
        REQUIRE_THAT(H1[2][3], WithinAbs(H[2][3], 1e-7));

        REQUIRE_THAT(H1[3][0], WithinAbs(H[3][0], 1e-7));
        REQUIRE_THAT(H1[3][1], WithinAbs(H[3][1], 1e-7));
        REQUIRE_THAT(H1[3][2], WithinAbs(H[3][2], 1e-7));
        REQUIRE_THAT(H1[3][3], WithinAbs(H[3][3], 1e-7));
    }

    SECTION("interpolation matrix is orthogonal")
    {
        constexpr auto N = 4;
        float t          = GENERATE(take(4, random(0.f, 1.f)));
        auto H           = dsp::hadamardInterpolMatrix<N>(t);

        // H is ortogonal
        REQUIRE_THAT(dot(H[0], H[0]), WithinAbs(1.f, 1e-6));
        REQUIRE_THAT(dot(H[0], H[1]), WithinAbs(0.f, 1e-6));
        REQUIRE_THAT(dot(H[0], H[2]), WithinAbs(0.f, 1e-6));
        REQUIRE_THAT(dot(H[0], H[3]), WithinAbs(0.f, 1e-6));

        REQUIRE_THAT(dot(H[1], H[1]), WithinAbs(1.f, 1e-6));
        REQUIRE_THAT(dot(H[1], H[2]), WithinAbs(0.f, 1e-6));
        REQUIRE_THAT(dot(H[1], H[3]), WithinAbs(0.f, 1e-6));

        REQUIRE_THAT(dot(H[2], H[2]), WithinAbs(1.f, 1e-6));
        REQUIRE_THAT(dot(H[2], H[3]), WithinAbs(0.f, 1e-6));

        REQUIRE_THAT(dot(H[3], H[3]), WithinAbs(1.f, 1e-6));
    }
}
