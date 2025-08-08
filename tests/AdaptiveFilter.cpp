#include "dsp/AdaptiveFilter.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch::Matchers;

template <template <typename, size_t> class Algo, typename T, size_t N = 500>
static void testSine()
{
    // as seen in
    // https://ocw.mit.edu/courses/2-161-signal-processing-continuous-and-discrete-fall-2008/resources/rls/
    dsp::AdaptiveFilter<T, 2> afilter;
    Algo<T, 2> rls{0.99};
    typename decltype(afilter)::DL delay;

    constexpr size_t kN = N;
    const auto &a       = afilter.getCoeffs();
    for (size_t n = 0; n < kN; ++n) {
        auto sig = std::sin(dsp::constants<T>::pi * 2 * T(n) / 12.f);
        rls.process(dsp::Context(&sig), delay, afilter);
    }
    REQUIRE_THAT(a.get(0), WithinAbs(-1.f, 1e-4f));
    REQUIRE_THAT(a.get(1), WithinAbs(std::sqrt(3.f), 1e-4f));
}

template <template <typename, size_t> class Algo, typename T, size_t N = 100>
static void testReconstruct()
{
    dsp::AdaptiveFilter<T, 2> afilter;
    Algo<T, 2> rls{0.99};
    typename decltype(afilter)::DL delay;
    typename decltype(afilter)::DL delayout;

    constexpr size_t kN = N;
    for (size_t n = 0; n < kN; ++n) {
        auto sig = T(rand()) / T(INT_MAX);
        sig += std::sin(dsp::constants<T>::pi * 2 * T(n) / 12.f);
        auto x = sig;
        dsp::Context ctxt{&sig};
        auto oafilter = afilter;
        rls.process(ctxt, delay, afilter);
        oafilter.reconstruct(ctxt, delayout);
        // reconstructed is the same as original
        REQUIRE_THAT(sig, WithinAbs(x, 1e-6f));
    }
}

template <typename T> static void testWIIRFilter0()
{
    // test if WIIR filter with warp arg = 0 is the same as normal IIR Filter
    T adata[] = {-0.5, 0.5};
    dsp::linalg::Vector<T, 1> a(adata);
    dsp::WarpedIIR wiir(a);

    dsp::CopyDelayLine<T, 1> dl{};

    typename decltype(wiir)::State wiirState{};
    T warp = 0.0;
    wiir.setCoeff(warp);
    dsp::AdaptiveFilter<T, 1> &iir = wiir;

    for (int i = 0; i < 30; ++i) {
        T x    = i == 0 ? 1.f : 0.f;
        T sig  = x;
        T sig2 = x;

        wiir.reconstruct(dsp::Context(&sig), wiirState);
        iir.reconstruct(dsp::Context(&sig2), dl);
        REQUIRE_THAT(sig, WithinAbs(sig2, 1e-6f));
    }
}

template <template <typename, size_t> class Algo, typename T, size_t N = 100>
static void testReconstructWarped()
{
    dsp::WarpedIIR<T, 2> afilter;
    Algo<T, 2> rls{0.99};
    typename decltype(afilter)::DL delay;
    typename decltype(afilter)::State rstate{};

    constexpr size_t kN = N;
    for (size_t n = 0; n < kN; ++n) {
        auto sig = T(rand()) / T(INT_MAX);
        sig += std::sin(dsp::constants<T>::pi * 2 * T(n) / 12.f);
        auto x = sig;
        dsp::Context ctxt{&sig};
        auto oafilter = afilter;
        rls.process(ctxt, delay, afilter);
        oafilter.reconstruct(ctxt, rstate);
        // reconstructed is the same as original
        REQUIRE_THAT(sig, WithinAbs(x, 1e-6f));
    }
}

TEST_CASE("Adaptive Filter", "[dsp][adaptive][filter]")
{
    SECTION("RLS")
    {
        SECTION("Sine")
        {
            SECTION("float") { testSine<dsp::RLS, float>(); }
            SECTION("double") { testSine<dsp::RLS, double>(); }
        }
        SECTION("Reconstruct")
        {
            SECTION("float") { testReconstruct<dsp::RLS, float>(); }
            SECTION("double") { testReconstruct<dsp::RLS, double>(); }
        }
    }
    SECTION("RLS-DCD")
    {
        SECTION("Sine")
        {
            SECTION("float") { testSine<dsp::RLSDCD, float>(); }
            SECTION("double") { testSine<dsp::RLSDCD, double>(); }
        }
        SECTION("Reconstruct")
        {
            SECTION("float") { testReconstruct<dsp::RLSDCD, float>(); }
            SECTION("double") { testReconstruct<dsp::RLSDCD, double>(); }
        }
    }

    SECTION("Warped IIR")
    {
        SECTION("float") { testWIIRFilter0<float>(); }
        SECTION("double") { testWIIRFilter0<double>(); }
    }

    SECTION("Reconstruct Warped")
    {
        SECTION("float") { testReconstructWarped<dsp::RLSDCD, float>(); }
        SECTION("double") { testReconstructWarped<dsp::RLSDCD, double>(); }
    }
}
