#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "dsp/FastMath.h"

using namespace Catch::Matchers;

template <typename T> static void testTanh()
{
    T eps      = std::is_same_v<T, double> ? 1e-6 : 1e-5;
    T values[] = {0, 0.1, 0.43, 0.83, 1, 1.238, 3.398, 5.4983, 7.3873};
    for (auto v : values) {
        REQUIRE_THAT(dsp::fasttanh(v), WithinAbs(std::tanh(v), eps));
        REQUIRE_THAT(dsp::fasttanh(-v), WithinAbs(std::tanh(-v), eps));
    }
}

TEST_CASE("Fast tanh", "[dsp][math][tanh]")
{
    SECTION("float") { testTanh<float>(); }
    SECTION("double") { testTanh<double>(); }
}
