#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "dsp/FastMath.h"

using namespace Catch::Matchers;

template <typename T> static void testTanh()
{
    T eps      = std::is_same_v<T, double> ? T(1e-6) : T(1e-5);
    T values[] = {T(0),     T(0.1),   T(0.43),   T(0.83),  T(1),
                  T(1.238), T(3.398), T(5.4983), T(7.3873)};
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
