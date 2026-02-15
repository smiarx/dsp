#include "dsp/FIRFilter.h"
#include "dsp/Buffer.h"
#include "dsp/Context.h"
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdlib>

template <typename T, size_t Order, bool Vec = false> static void testFir()
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

template <typename T, size_t Order, size_t M, bool Vec = false>
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
            auto val            = dsp::get(xDecimated[n / M], i);
            constexpr auto kEps = decltype(val)(1e-6);
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

    SECTION("python data")
    {
        constexpr auto kN       = 50;
        constexpr auto kO       = 13;
        std::array<float, kN> x = {
            0.f,         0.01923305f, 0.03845899f, 0.05767071f, 0.07686108f,
            0.09602303f, 0.11514945f, 0.13423327f, 0.15326743f, 0.17224489f,
            0.19115863f, 0.21000165f, 0.22876698f, 0.24744769f, 0.26603685f,
            0.28452759f, 0.30291307f, 0.32118649f, 0.33934109f, 0.35737015f,
            0.375267f,   0.39302503f, 0.41063766f, 0.42809838f, 0.44540072f,
            0.46253829f, 0.47950475f, 0.49629381f, 0.51289928f, 0.529315f,
            0.5455349f,  0.56155299f, 0.57736333f, 0.59296008f, 0.60833747f,
            0.6234898f,  0.63841148f, 0.65309698f, 0.66754088f, 0.68173782f,
            0.69568255f, 0.70936992f, 0.72279486f, 0.73595241f, 0.7488377f,
            0.76144596f, 0.77377252f, 0.78581284f, 0.79756244f, 0.80901699f,
        };

        // signal.decimate(x, 2,27,'fir', zero_phase=False)
        std::array<float, kN / 2> expect = {
            0.00000000e+00f,  8.30978079e-05f,  2.40111901e-05f,
            3.10682470e-04f,  -1.27662335e-04f, 8.51615654e-04f,
            -1.13936683e-03f, 9.48246747e-03f,  4.87176327e-02f,
            8.64312536e-02f,  1.24966445e-01f,  1.62826426e-01f,
            2.00664657e-01f,  2.38123860e-01f,  2.75271734e-01f,
            3.12037896e-01f,  3.48342354e-01f,  3.84131390e-01f,
            4.19352050e-01f,  4.53952219e-01f,  4.87880701e-01f,
            5.21087295e-01f,  5.53522866e-01f,  5.85139421e-01f,
            6.15890181e-01f,
        };

        dsp::FIRDecimate<float, kO, 2, dsp::windows::Hamming> decimate;
        decltype(decimate)::DL decState;

        // delayline buffer
        dsp::Buffer<float, nextTo(decState)> buffer{};
        std::array<float, decltype(buffer)::kSize> bufdata{};
        buffer.setData(bufdata.data());

        std::array<float, kN> xDec{};

        dsp::BufferContext ctxtIn(x.data(), kN, buffer);
        dsp::Context ctxtOut(xDec.data(), kN);

        decimate.decimate(ctxtIn, ctxtOut, decState, 0);

        for (size_t i = 0; i < size_t(ctxtOut.getBlockSize()); ++i) {
            REQUIRE_THAT(xDec[i], Catch::Matchers::WithinAbs(expect[i], 1e-6f));
        }
    }
}

template <typename T, size_t Order, size_t L, bool Vec = false>
static void testFirInterpolate()
{
    // build decimate filter
    dsp::FIRInterpolate<T, Order, L> interpolate{};
    typename decltype(interpolate)::DL interpState;

    // reproduce interpolate filter coeffs
    constexpr size_t kNCoeffs = (Order + 1) * L;
    std::array<T, kNCoeffs> bF;
    auto freq = dsp::baseType<T>(1) / L;
    auto mid  = (kNCoeffs - 1) / dsp::baseType<T>(2);
    for (size_t i = 0; i < kNCoeffs; ++i) {
        auto fi = dsp::baseType<T>(i);
        bF[i]   = dsp::windows::Kaiser<140>::generate((fi - mid) / mid) *
                dsp::sinc((fi - mid) * freq);
    }
    // scale
    auto sum = dsp::load(T(0));
    for (auto &b : bF) sum += b;
    auto scale = T(L) / sum;
    for (auto &b : bF) b = b * scale;

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
            auto val            = dsp::get(xInterp[n], i);
            constexpr auto kEps = decltype(val)(1e-6);
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
