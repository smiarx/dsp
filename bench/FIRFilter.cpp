#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include "dsp/Buffer.h"
#include "dsp/Context.h"
#include "dsp/FIRFilter.h"

TEST_CASE("FIR Filter")
{
    constexpr size_t kK = 2;
    constexpr size_t kN = 512;

    constexpr size_t kOrder = 15;
    dsp::fSample<kK> x[kN];
    dsp::fSample<kK> x2[kN];

    for (auto &xn : x)
        for (size_t i = 0; i < kK; ++i)
            xn[i] = GENERATE(take(1, random(-1.f, 1.f)));
    for (auto &x2n : x2)
        for (size_t i = 0; i < kK; ++i)
            x2n[i] = GENERATE(take(1, random(-1.f, 1.f)));

    std::array<dsp::fData<kK>, kOrder + 1> b;
    for (auto &b0 : b)
        for (size_t i = 0; i < kK; ++i)
            b0[i] = GENERATE(take(1, random(-1.f, 1.f)));

    auto ctxt = dsp::Context(x, kN);
    dsp::FIRFilter<kK, kOrder> filter(b);
    decltype(filter)::DL filterState{};

    BENCHMARK("FIR filter")
    {
        contextFor(ctxt) { filter.process(c, filterState); }
        return x;
    };

    dsp::FIRDecimate<kK, kOrder, 3> decimate;
    dsp::FIRInterpolate<kK, kOrder, 3> interpolate;

    decltype(decimate)::DL<0> dldecimate;
    decltype(interpolate)::DL dlinterpolate;

    dsp::Buffer<dsp::fSample<kK>, nextTo(dldecimate)> buffer;
    dsp::fSample<kK> bufdata[decltype(buffer)::kSize];
    buffer.setBuffer(bufdata);

    auto bufctxt    = dsp::BufferContext(x, kN, buffer);
    auto ctxtDec    = dsp::Context(x2, kN);
    auto ctxtInterp = dsp::Context(x, kN);

    BENCHMARK("FIR decimate")
    {
        decimate.decimate(bufctxt, ctxtDec, dldecimate, 0);
        return x;
    };
    BENCHMARK("FIR interpolate")
    {
        interpolate.interpolate(ctxtDec, ctxtInterp, dlinterpolate, 0);
        return x;
    };
}
