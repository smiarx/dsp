#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include "dsp/Context.h"
#include "dsp/FIRFilter.h"

TEST_CASE("FIR Filter")
{
    constexpr size_t kN     = 2;
    constexpr size_t kOrder = 26;
    dsp::fSample<kN> x;
    for (size_t i = 0; i < kN; ++i) {
        x[i] = GENERATE(take(1, random(-1.f, 1.f)));
    }

    std::array<dsp::fData<kN>, kOrder + 1> b;
    for (auto &b0 : b)
        for (size_t i = 0; i < kN; ++i)
            b0[i] = GENERATE(take(1, random(-1.f, 1.f)));

    auto ctxt = dsp::Context(&x, 1);

    dsp::FIRFilter<kN, kOrder> filter(b);
    decltype(filter)::DL filterState{};

    BENCHMARK("FIR Order 11 process")
    {
        filter.process(ctxt, filterState);
        return x;
    };
}
