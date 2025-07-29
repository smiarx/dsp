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
    using ft            = dsp::mfloat<kK>;

    ft x[kN];
    for (auto &xn : x)
        for (size_t i = 0; i < kK; ++i) {
            xn[i] = GENERATE(take(1, random(-1.f, 1.f)));
        }

    float a0 = GENERATE(take(1, random(-1.f, 1.f)));
    float a1 = GENERATE(take(1, random(-1.f, 1.f)));

    auto ctxt = dsp::Context(x, kN);

    {
        dsp::AllPass<ft> ap1;
        dsp::CopyDelayLine<ft, 1> ap1state{};
        ap1.setCoeff({a0, a0});
        BENCHMARK("1st order")
        {
            CTXTRUN(ctxt) { ap1.process(ctxt, ap1state); };
            return x;
        };
    }
    {
        dsp::AllPass2<ft> ap2;
        decltype(ap2)::State ap2state{};
        ap2.setCoeffs({a0, a0}, {a1, a1});
        BENCHMARK("2nd order")
        {
            CTXTRUN(ctxt) { ap2.process(ctxt, ap2state); };
            return x;
        };
    }
}
