#include "dsp/FIRFilter.h"
#include <catch2/catch_test_macros.hpp>

template <int N, int Order, bool Vec = false> void test_fir()
{
    std::array<dsp::Data<float, N>, Order + 1> b;
    arrayFor(b, k) { arrayFor(b[0], i) b[k][i] = rand(); }
    dsp::FIRFilter<N, Order> filter(b);
    typename decltype(filter)::DL fstate;

    std::array<dsp::Sample<float, N>, (Order + 1) + 8> x = {0.f};
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
    test_fir<1, 8>();
    test_fir<1, 47>();
    test_fir<2, 23>();
    test_fir<4, 17>();
    test_fir<1, 8, true>();
    test_fir<1, 47, true>();
    test_fir<2, 23, true>();
    test_fir<4, 15, true>();
}
