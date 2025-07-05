#include "dsp/FIRFilter.h"
#include "dsp/Context.h"
#include "dsp/Signal.h"
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstdlib>

template <typename T, int Order, bool Vec = false> static void testFir()
{
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

    std::array<T, (Order + 1) + 8> x = {};
    x[0]                             = 1;
    dsp::Context<T, Vec> ctxt(x.data(), Order + 2);

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
