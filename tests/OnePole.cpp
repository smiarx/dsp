#include "dsp/OnePole.h"
#include "dsp/Context.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("OnePole", "[dsp][onepole]")
{
    SECTION("scalar")
    {
        int N       = GENERATE(100, 50, 300);
        double init = GENERATE(take(3, random(-10., 10.)));
        double goal = GENERATE(take(3, random(-10., 10.)));

        constexpr auto kConvergenceDiff = 0.01f;
        dsp::OnePole<double> ev;
        ev.setRate(kConvergenceDiff, (1.f / static_cast<double>(N)));
        decltype(ev)::State evstate = init;

        double x;
        for (int n = 0; n < N; ++n) {
            x = goal;
            ev.process(dsp::Context(&x), evstate);
        }
        REQUIRE_THAT(x, Catch::Matchers::WithinAbs(
                            goal - (goal - init) * kConvergenceDiff, 1e-6));
    }

    SECTION("vector")
    {
        using ft = dsp::mfloat<4>;

        dsp::OnePole<ft> ev;
        ft convdiff = {0.01,0.08,0.001,0.123};
        int N[] = {100,90,50,74};
        ft rate;
        for(size_t i = 0; i < 4; ++i)
            rate[i] = 1.f/static_cast<float>(N[i]);

        ev.setRate(convdiff, rate);
        decltype(ev)::State evstate = {};

        ft x;
        for (int n = 0; n < 101; ++n) {
            for(int i = 0; i < 4; ++i)
            {
                if(n == N[i])
                {
                    REQUIRE_THAT(x[i], Catch::Matchers::WithinAbs(
                                        1 - convdiff[i], 1e-6));
                }
            }
            x = {1,1,1,1};
            ev.process(dsp::Context(&x), evstate);
        }
    }
}
