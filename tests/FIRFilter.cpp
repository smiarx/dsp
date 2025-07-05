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
            REQUIRE_THAT(
                dsp::get(xDecimated[n / M], i),
                Catch::Matchers::WithinRel(dsp::get(xFilt[n], i), 1e-6f));
        }
    }
}

TEST_CASE("FIR decimate", "[dsp][fir][decimate]")
{
    SECTION("float x 1") { testFirDecimate<float, 15, 3>(); }
    SECTION("float x 2") { testFirDecimate<dsp::mfloat<2>, 19, 2>(); }
    SECTION("double x 2") { testFirDecimate<dsp::mfloat<2>, 12, 4>(); }
    SECTION("float x 4") { testFirDecimate<dsp::mfloat<4>, 23, 5>(); }
    SECTION("float x 1 vec") { testFirDecimate<float, 15, 3, true>(); }
    SECTION("float x 2 vec") { testFirDecimate<dsp::mfloat<2>, 19, 2, true>(); }
    SECTION("float x 4 vec") { testFirDecimate<dsp::mfloat<4>, 23, 5, true>(); }
}
