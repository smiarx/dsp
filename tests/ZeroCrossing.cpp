#include "dsp/ZeroCrossing.h"
#include "dsp/Context.h"
#include <catch2/catch_all.hpp>

TEST_CASE("Zero Crossing", "[dsp][zerocrossing]")
{
    SECTION("scalar")
    {
        dsp::ZeroCrossing<float> zc;

        float values[] = {0, 1, 2, 3, -1, -4, 4, -2, -4};
        int expect[]   = {0, 0, 0, 0, 1, 0, 1, 1, 0};

        dsp::Context ctxtv(values, sizeof(values) / sizeof(values[0]));
        dsp::Context ctxte(expect, sizeof(expect) / sizeof(expect[0]));

        CTXTRUN2(ctxtv, ctxte)
        {
            REQUIRE(zc.process(ctxtv) == ctxte.getInput());
        };
    }

    SECTION("vector")
    {
        using ft = dsp::mfloat<2>;
        dsp::ZeroCrossing<ft> zc;

        ft values[] = {{0, -1}, {-2, -2}, {-2, 3}, {-1, 3}, {-4, 4}, {-2, -4}};
        dsp::mint<2> expect[] = {{0, 1}, {1, 0}, {0, 1},
                                 {0, 0}, {0, 0}, {0, 1}};

        dsp::Context ctxtv(values, sizeof(values) / sizeof(values[0]));
        dsp::Context ctxte(expect, sizeof(expect) / sizeof(expect[0]));

        CTXTRUN2(ctxtv, ctxte)
        {
            auto r = zc.process(ctxtv);
            auto e = ctxte.getInput();
            REQUIRE(r[0] == e[0]);
            REQUIRE(r[1] == e[1]);
        };
    }
}
