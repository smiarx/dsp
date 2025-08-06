#include "dsp/AdaptiveFilter.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch::Matchers;

template <template <typename, size_t> class Algo, typename T, size_t N = 500>
static void testSine()
{
    // as seen in
    // https://ocw.mit.edu/courses/2-161-signal-processing-continuous-and-discrete-fall-2008/resources/rls/
    dsp::AdaptiveFilter<T, 2> afilter;
    Algo<T, 2> rls{0.99};
    typename decltype(afilter)::DL delay;

    constexpr size_t kN = N;
    const auto &a       = afilter.getCoeffs();
    for (size_t n = 0; n < kN; ++n) {
        auto sig = std::sin(dsp::constants<T>::pi * 2 * T(n) / 12.f);
        rls.process(dsp::Context(&sig), delay, afilter);
    }
    REQUIRE_THAT(a.get(0), WithinAbs(-1.f, 1e-4f));
    REQUIRE_THAT(a.get(1), WithinAbs(std::sqrt(3.f), 1e-4f));
}

TEST_CASE("Adaptive Filter", "[dsp][adaptive][filter]")
{
    SECTION("RLS")
    {
        SECTION("Sine")
        {
            SECTION("float") { testSine<dsp::RLS, float>(); }
            SECTION("double") { testSine<dsp::RLS, double>(); }
        }
    }
    SECTION("RLS-DCD")
    {
        SECTION("Sine")
        {
            SECTION("float") { testSine<dsp::RLSDCD, float>(); }
            SECTION("double") { testSine<dsp::RLSDCD, double>(); }
        }
    }
}
