#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include "dsp/Context.h"
#include "dsp/VAFilters.h"

#define BENCHMARK_FILTER(Name, Type)                        \
    {                                                       \
        dsp::va::Name<kK, dsp::va::k##Type> filter;         \
        decltype(filter)::State state;                      \
        auto freq = GENERATE(take(1, random(0.f, 1.f)));    \
        filter.setFreq({freq, freq});                       \
        BENCHMARK(#Type)                                    \
        {                                                   \
            contextFor(ctxt) { filter.process(c, state); }; \
        };                                                  \
    }

constexpr size_t kK     = 2;
constexpr size_t kN     = 512;
constexpr size_t kOrder = 26;
using ft                = dsp::fSample<kK>;
static ft x[kN];
TEST_CASE("VA Filters", "[dsp][vafilters]")
{
    for (auto &xn : x)
        for (size_t i = 0; i < kK; ++i)
            xn[i] = GENERATE(take(1, random(-1.f, 1.f)));

    auto ctxt = dsp::Context(x, kN);

    SECTION("OnePole")
    {
        BENCHMARK_FILTER(OnePole, LowPass);
        BENCHMARK_FILTER(OnePole, HighPass);
        BENCHMARK_FILTER(OnePole, AllPass);
    }
    SECTION("SVF")
    {
        BENCHMARK_FILTER(SVF, LowPass);
        BENCHMARK_FILTER(SVF, HighPass);
        BENCHMARK_FILTER(SVF, AllPass);
        BENCHMARK_FILTER(SVF, BandPass);
        BENCHMARK_FILTER(SVF, Notch);
    }
    SECTION("Ladder")
    {
        BENCHMARK_FILTER(Ladder, LowPass);
        BENCHMARK_FILTER(Ladder, HighPass);
    }
}
