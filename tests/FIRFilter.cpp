#include "dsp/FIRFilter.h"
#include "dsp/Buffer.h"
#include "dsp/Context.h"
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdlib>

template <typename T, int Order, bool Vec = false> static void testFir()
{
    // generate random coeffs
    std::array<T, Order + 1> b;
    for (size_t k = 0; k < Order + 1; ++k)
        if constexpr (dsp::kTypeWidth<T> > 1) {
            for (size_t i = 0; i < dsp::kTypeWidth<T>; ++i)
                b[k][i] = static_cast<float>(std::rand()) /
                          static_cast<float>(INT32_MAX);
        } else
            b[k] =
                static_cast<float>(std::rand()) / static_cast<float>(INT32_MAX);

    dsp::FIRFilter<T, Order> filter(b);
    typename decltype(filter)::DL fstate;

    // input is impulse
    std::array<T, (Order + 1) + 8> x = {};
    x[0]                             = 1;
    dsp::Context<T, Vec> ctxt(x.data(), Order + 2);

    // check if output is impulse response given FIR coeffs
    size_t n = 0;
    CTXTRUN(ctxt)
    {
        filter.process(ctxt, fstate);
        auto xout = ctxt.getInput();

        for (size_t k = 0; k < decltype(ctxt)::kIncrSize; ++k) {
            for (size_t i = 0; i < dsp::kTypeWidth<T>; ++i) {
                auto ki = k * dsp::kTypeWidth<T> + i;
                if (n + k < Order + 1)
                    REQUIRE(dsp::get(xout, ki) == dsp::get(b[n + k], i));
                else {
                    REQUIRE(dsp::get(xout, ki) == 0.f);
                }
            }
        }
        n += decltype(ctxt)::kIncrSize;
    };
}

TEST_CASE("FIR filter test", "[dsp][firfilter]")
{
    SECTION("float x 1") { testFir<float, 8>(); }
    SECTION("double x 1") { testFir<double, 47>(); }
    SECTION("float x 2") { testFir<dsp::mfloat<2>, 25>(); }
    SECTION("double x 2") { testFir<dsp::mdouble<2>, 23>(); }
    SECTION("float x 4") { testFir<dsp::mfloat<4>, 17>(); }
    SECTION("float x 1 vec") { testFir<float, 8, true>(); }
    SECTION("double x 1 vec") { testFir<double, 47, true>(); }
    SECTION("float x 2 vec") { testFir<dsp::mfloat<2>, 23, true>(); }
    SECTION("float x 4 vec") { testFir<dsp::mfloat<4>, 15, true>(); }
}

template <typename T, int Order, int M, bool Vec = false>
static void testFirDecimate()
{
    // build decimate filter
    dsp::FIRDecimate<T, Order, M> decimate{};
    typename decltype(decimate)::DL decState;

    auto b                      = decimate.getCoeffs();
    constexpr size_t kFiltOrder = (Order + 1) * M - 1;
    constexpr auto kPad         = decltype(decimate)::kPad;

    // copy decimate coeffs to filter coeffs
    std::array<T, kFiltOrder + 1> bF;
    std::copy(std::begin(b) + kPad, std::end(b) - kPad + 1, std::begin(bF));
    std::reverse(std::begin(bF), std::end(bF));

    dsp::FIRFilter<T, kFiltOrder> filter(bF);
    typename decltype(filter)::DL filtState;

    // processed data
    constexpr size_t kNSamples          = 41;
    std::array<T, kNSamples> xDecimated = {};
    std::array<T, kNSamples> xOrig      = {};
    std::array<T, kNSamples> xFilt      = {};
    xOrig[0] = xFilt[0] = 1;
    xOrig[3] = xFilt[3] = 1;
    xOrig[4] = xFilt[4] = 1;

    // delayline buffer
    dsp::Buffer<T, nextTo(decState)> buffer{};
    std::array<T, decltype(buffer)::kSize> bufdata{};
    buffer.setData(bufdata.data());

    // run filter first
    dsp::Context<T, Vec> ctxtFilt(xFilt.data(), kNSamples);
    CTXTRUN(ctxtFilt) { filter.process(ctxtFilt, filtState); };

    // decimate
    dsp::BufferContext<T, nextTo(decState), Vec> ctxtIn(xOrig.data(), kNSamples,
                                                        buffer);
    dsp::Context<T> ctxtOut(xDecimated.data(), kNSamples);
    decimate.decimate(ctxtIn, ctxtOut, decState, 0);

    // check that decimate is equal to 1 over M filtered sample
    for (size_t n = 0; n < kNSamples; n += M) {
        for (size_t i = 0; i < dsp::kTypeWidth<T>; ++i) {
            auto val                     = dsp::get(xDecimated[n / M], i);
            constexpr decltype(val) kEps = 1e-6;
            REQUIRE_THAT(
                val, Catch::Matchers::WithinRel(dsp::get(xFilt[n], i), kEps));
        }
    }
}

TEST_CASE("FIR decimate", "[dsp][fir][decimate]")
{
    SECTION("float x 1") { testFirDecimate<float, 15, 3>(); }
    SECTION("float x 2") { testFirDecimate<dsp::mfloat<2>, 19, 2>(); }
    SECTION("double x 2") { testFirDecimate<dsp::mdouble<2>, 12, 4>(); }
    SECTION("float x 4") { testFirDecimate<dsp::mfloat<4>, 23, 5>(); }
    SECTION("float x 1 vec") { testFirDecimate<float, 15, 3, true>(); }
    SECTION("float x 2 vec") { testFirDecimate<dsp::mfloat<2>, 19, 2, true>(); }
    SECTION("float x 4 vec") { testFirDecimate<dsp::mfloat<4>, 23, 5, true>(); }
}

template <typename T, int Order, int L, bool Vec = false>
static void testFirInterpolate()
{
    // build decimate filter
    dsp::FIRInterpolate<T, Order, L> interpolate{};
    typename decltype(interpolate)::DL interpState;

    // reproduce interpolate filter coeffs
    constexpr auto kNCoeffs = (Order + 1) * L;
    std::array<T, kNCoeffs> bF;
    double freq = 1. / L;
    double mid  = (kNCoeffs - 1) / 2.f;
    for (auto i = 0; i < kNCoeffs; ++i) {
        double fi = i;
        bF[i]     = dsp::window::Kaiser<140>::generate((fi - mid) / mid) *
                dsp::sinc((fi - mid) * freq);
    }

    constexpr auto kFiltOrder = kNCoeffs - 1;
    dsp::FIRFilter<T, kFiltOrder> filter(bF);
    typename decltype(filter)::DL filtState;

    // processed data
    constexpr size_t kNSamples           = 41;
    std::array<T, kNSamples * L> xInterp = {};
    std::array<T, kNSamples> xOrig       = {};
    std::array<T, kNSamples * L> xFilt   = {};
    xOrig[0] = xFilt[0] = 1;
    xOrig[3] = xFilt[3 * L] = 1;
    xOrig[4] = xFilt[4 * L] = 1;

    // delayline buffer
    dsp::Buffer<T, nextTo(interpState)> buffer{};
    std::array<T, decltype(buffer)::kSize> bufdata{};
    buffer.setData(bufdata.data());

    // run filter first
    dsp::Context<T, Vec> ctxtFilt(xFilt.data(), kNSamples * L);
    CTXTRUN(ctxtFilt) { filter.process(ctxtFilt, filtState); };

    // interpolate
    dsp::BufferContext<T, nextTo(interpState), Vec> ctxtIn(xOrig.data(),
                                                           kNSamples, buffer);
    // dsp::Context<T, Vec> ctxtIn(xOrig.data(), kNSamples);
    dsp::Context<T> ctxtOut(xInterp.data(), kNSamples * L);
    interpolate.interpolate(ctxtIn, ctxtOut, interpState, 0);

    // check that decimate is equal to 1 over M filtered sample
    for (size_t n = 0; n < kNSamples * L; ++n) {
        for (size_t i = 0; i < dsp::kTypeWidth<T>; ++i) {
            auto val                     = dsp::get(xInterp[n], i);
            constexpr decltype(val) kEps = 1e-6;
            REQUIRE_THAT(
                val, Catch::Matchers::WithinAbs(dsp::get(xFilt[n], i), kEps));
        }
    }
}

TEST_CASE("FIR interpolate", "[dsp][fir][interpolate]")
{
    SECTION("float x 1") { testFirInterpolate<float, 15, 3>(); }
    SECTION("float x 2") { testFirInterpolate<dsp::mfloat<2>, 19, 2>(); }
    SECTION("double x 2") { testFirInterpolate<dsp::mdouble<2>, 12, 4>(); }
    SECTION("float x 4") { testFirInterpolate<dsp::mfloat<4>, 23, 5>(); }
}
