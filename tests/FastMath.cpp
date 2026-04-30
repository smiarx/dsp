#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dsp/math/FastMath.hpp>

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

template <typename T> static void testExp2()
{
    constexpr T kEps = T(1e-4);
    T values[]       = {T(0),     T(0.1),   T(0.43),   T(0.83),  T(1),
                        T(1.238), T(3.398), T(5.4983), T(7.3873)};
    for (auto v : values) {
        REQUIRE_THAT(dsp::fastexp2(v), WithinRel(std::exp2(v), kEps));
        REQUIRE_THAT(dsp::fastexp2(-v), WithinRel(std::exp2(-v), kEps));
    }
}
TEST_CASE("Fast exp2", "[dsp][math][exp2]")
{
    SECTION("float") { testExp2<float>(); }
    SECTION("double") { testExp2<double>(); }
    SECTION("simd float")
    {
        dsp::mfloat4 x = {0.1f, 1.f, -0.4f, 3.f};
        auto r         = dsp::fastexp2(dsp::load(x));

        for (size_t i = 0; i < 4; ++i) {
            REQUIRE_THAT(r[i], WithinRel(std::exp2(x[i]), 1e-4f));
        }
    }
}

template <typename T> static void testLog2()
{
    constexpr T kEps = T(2e-4);
    T values[]       = {T(0.1),   T(0.43),  T(0.83),   T(1),
                        T(1.238), T(3.398), T(5.4983), T(7.3873)};
    for (auto v : values) {
        REQUIRE_THAT(dsp::fastlog2(v), WithinAbs(std::log2(v), kEps));
    }
}
TEST_CASE("Fast log2", "[dsp][math][log2]")
{
    SECTION("float") { testLog2<float>(); }
    SECTION("double") { testLog2<double>(); }
    SECTION("simd float")
    {
        dsp::mfloat4 x = {0.1f, 1.f, 3.4f, 8.2f};
        auto r         = dsp::fastlog2(dsp::load(x));

        for (size_t i = 0; i < 4; ++i) {
            REQUIRE_THAT(r[i], WithinAbs(std::log2(x[i]), 2e-4f));
        }
    }
}
