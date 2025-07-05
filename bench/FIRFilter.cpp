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
    constexpr size_t kK     = 2;
    constexpr size_t kN     = 512;
    using ft                = dsp::mfloat<kK>;
    constexpr size_t kOrder = 15;
    ft x[kN];
    ft x2[kN];

    for (auto &xn : x)
        for (size_t i = 0; i < kK; ++i)
            xn[i] = GENERATE(take(1, random(-1.f, 1.f)));
    for (auto &x2n : x2)
        for (size_t i = 0; i < kK; ++i)
            x2n[i] = GENERATE(take(1, random(-1.f, 1.f)));

    std::array<ft, kOrder + 1> b;
    for (auto &b0 : b)
        for (size_t i = 0; i < kK; ++i)
            b0[i] = GENERATE(take(1, random(-1.f, 1.f)));

    auto ctxt = dsp::Context(x, kN);
    dsp::FIRFilter<ft, kOrder> filter(b);
    decltype(filter)::DL filterState{};

    BENCHMARK("FIR filter")
    {
        CTXTRUN(ctxt) { filter.process(ctxt, filterState); };
        return x;
    };

    /*
    dsp::FIRDecimate<ft, kOrder, 3> decimate;
    dsp::FIRInterpolate<ft, kOrder, 3> interpolate;

    decltype(decimate)::DL<0> dldecimate;
    decltype(interpolate)::DL<nextTo(dldecimate)> dlinterpolate;

    dsp::Buffer<dsp::fSample<kK>, nextTo(dlinterpolate)> buffer;
    dsp::fSample<kK> bufdata[decltype(buffer)::kSize];
    buffer.setBuffer(bufdata);

    auto bufctxt    = dsp::BufferContext(x, kN, buffer);
    auto ctxtDec    = dsp::BufferContext(x2, kN, buffer);
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
    */
}
