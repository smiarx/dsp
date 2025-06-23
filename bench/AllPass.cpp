#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include "dsp/AllPass.h"
#include "dsp/Context.h"

TEST_CASE("Allpass")
{
    constexpr size_t kK = 2;
    constexpr size_t kN = 512;

    dsp::fSample<kK> x[kN];
    for (auto &xn : x)
        for (size_t i = 0; i < kK; ++i) {
            xn[i] = GENERATE(take(1, random(-1.f, 1.f)));
        }

    float a0 = GENERATE(take(1, random(-1.f, 1.f)));
    float a1 = GENERATE(take(1, random(-1.f, 1.f)));

    auto ctxt = dsp::Context(x, kN);

    {
        dsp::AllPass<kK> ap1;
        dsp::CopyDelayLine<kK, 1> ap1state{};
        ap1.setCoeff({a0, a0});
        BENCHMARK("1st order")
        {
            contextFor(ctxt) { ap1.process(c, ap1state); }
            return x;
        };
    }
    {
        dsp::AllPass2<kK> ap2;
        decltype(ap2)::State ap2state{};
        ap2.setCoeffs({a0, a0}, {a1, a1});
        BENCHMARK("2nd order")
        {
            contextFor(ctxt) { ap2.process(c, ap2state); }
            return x;
        };
    }
}
