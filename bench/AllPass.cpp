#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include "dsp/AllPass.h"
#include "dsp/Context.h"

TEST_CASE("Allpass")
{
    constexpr size_t kN = 2;
    dsp::fSample<kN> x;
    for (size_t i = 0; i < kN; ++i) {
        x[i] = GENERATE(take(1, random(-1.f, 1.f)));
    }

    float a0 = GENERATE(take(1, random(-1.f, 1.f)));
    float a1 = GENERATE(take(1, random(-1.f, 1.f)));

    auto ctxt = dsp::Context(&x, 1);

    {
        dsp::AllPass<kN> ap1;
        dsp::CopyDelayLine<kN, 1> ap1state{};
        ap1.setCoeff({a0, a0});
        BENCHMARK("1st order")
        {
            ap1.process(ctxt, ap1state);
            return x;
        };
    }
    {
        dsp::AllPass2<kN> ap2;
        decltype(ap2)::State ap2state{};
        ap2.setCoeffs({a0, a0}, {a1, a1});
        BENCHMARK("2nd order")
        {
            ap2.process(ctxt, ap2state);
            return x;
        };
    }
}
