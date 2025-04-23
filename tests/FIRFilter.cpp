#include "dsp/FIRFilter.h"
#include "dsp/Context.h"
#include "dsp/Signal.h"
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstdlib>

template <int N, int Order, bool Vec = false> static void testFir()
{
    std::array<dsp::Data<float, N>, Order + 1> b;
    arrayFor(b, k) { arrayFor(b[0], i) b[k][i] = rand(); }
    dsp::FIRFilter<N, Order> filter(b);
    typename decltype(filter)::DL fstate;

    std::array<dsp::Sample<float, N>, (Order + 1) + 8> x = {};
    arrayFor(x[0], i) { x[0][i] = 1.f; }
    dsp::Context<dsp::Sample<float, N>, Vec> ctxt(x.data(), Order + 2);

    contextFor(ctxt)
    {
        auto &xin = c.getSignal();
        filter.process(c, fstate);
        inFor(xin, k, i)
        {
            if (n + k < Order + 1) REQUIRE(xin[k][i] == b[n + k][i]);
            else {
                REQUIRE(xin[k][i] == 0.f);
            }
        }
    }
}

TEST_CASE("FIR filter test", "[dsp][firfilter]")
{
    testFir<1, 8>();
    testFir<1, 47>();
    testFir<2, 23>();
    testFir<4, 17>();
    testFir<1, 8, true>();
    testFir<1, 47, true>();
    testFir<2, 23, true>();
    testFir<4, 15, true>();
}
