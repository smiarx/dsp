#include "dsp/Context.h"
#include "dsp/Signal.h"
#include <array>
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Context Class", "[context]")
{
    SECTION("Construction and Copying")
    {
        dsp::Sample<float, 4> data;
        dsp::Context<dsp::Sample<float, 4>> ctx(&data);
        REQUIRE(ctx.getStep() == 1);
        REQUIRE(ctx.vecSize() == 1);

        // Copy construction
        dsp::Context<dsp::Sample<float, 4>> ctx_copy(ctx);
        REQUIRE(ctx_copy.getStep() == 1);
        REQUIRE(ctx_copy.vecSize() == 1);
    }

    SECTION("Vectorized Mode")
    {
        constexpr auto N = 2;
        dsp::Sample<float, N> data;
        dsp::Context<dsp::Sample<float, N>, true> ctx(&data);
        REQUIRE(ctx.getStep() == SIMDSIZE / sizeof(float) / N);
        REQUIRE(ctx.vecSize() == SIMDSIZE / sizeof(float) / N);
    }

    SECTION("Block Size and Increment Operations")
    {
        dsp::Sample<float, 4> data;
        dsp::Context<dsp::Sample<float, 4>> ctx(&data, 10);
        REQUIRE(ctx.getBlockSize() == 10);
        ctx.setBlockSize(5);
        REQUIRE(ctx.getBlockSize() == 5);

        ctx.next();
        REQUIRE(ctx.getBlockSize() == 5);
        ctx.nextBlock();
        REQUIRE(ctx.getBlockSize() == 5);
    }

    SECTION("Next and NextBlock Methods")
    {
        dsp::Sample<float, 4> data;
        dsp::Context<dsp::Sample<float, 4>> ctx(&data, 10);
        REQUIRE(ctx.getBlockSize() == 10);
        ctx.next(5);
        REQUIRE(ctx.getBlockSize() == 10);
        ctx.nextBlock();
        REQUIRE(ctx.getBlockSize() == 10);
    }
}

TEST_CASE("Helper Functions", "[helpers]")
{
    SECTION("_ctxtInfos")
    {
        dsp::Sample<float, 4> data1;
        dsp::Sample<float, 4> data2;
        dsp::Context<dsp::Sample<float, 4>> ctx1(&data1, 10);
        dsp::Context<dsp::Sample<float, 4>> ctx2(&data2, 10);

        int blockSize, incr;
        dsp::_ctxtInfos(blockSize, incr, ctx1, ctx2);
        REQUIRE(blockSize == 10);
        REQUIRE(incr == 1);
    }

    SECTION("_processBlock")
    {
        dsp::Sample<float, 4> data1;
        dsp::Sample<float, 4> data2;
        dsp::Context<dsp::Sample<float, 4>> ctx1(&data1, 10);
        dsp::Context<dsp::Sample<float, 4>> ctx2(&data2, 10);

        int count    = 0;
        auto process = [&](auto &c1, auto &c2) {
            count++;
        };

        dsp::_processBlock(process, ctx1, ctx2);
        REQUIRE(count == 10);
    }
}

TEST_CASE("Macros", "[macros]")
{
    SECTION("contextFor")
    {
        dsp::Sample<float, 4> data[10] = {
            {1.f, 2.f, 3.f, 4.f},     {5.f, 6.f, 7.f, 8.f},
            {9.f, 10.f, 11.f, 12.f},  {13.f, 14.f, 15.f, 16.f},
            {17.f, 18.f, 19.f, 20.f}, {21.f, 22.f, 23.f, 24.f},
            {25.f, 26.f, 27.f, 28.f}, {29.f, 30.f, 31.f, 32.f},
            {33.f, 34.f, 35.f, 36.f}, {37.f, 38.f, 39.f, 40.f},
        };
        dsp::Context<dsp::Sample<float, 4>> ctx(data, 10);
        int count = 0;
        contextFor(ctx)
        {
            if (count == 4) {
                auto &x = c.getSignal();
                REQUIRE(x[0][0] == 17.f);
                REQUIRE(x[0][1] == 18.f);
                REQUIRE(x[0][2] == 19.f);
                REQUIRE(x[0][3] == 20.f);
            }
            count++;
        }
        REQUIRE(count == 10);
    }

    SECTION("vecFor")
    {
        dsp::Sample<float, 4> data;
        dsp::Context<dsp::Sample<float, 4>> ctx(&data);
        int count = 0;
        vecFor(ctx) { count++; }
        REQUIRE(count == 1);
    }

    SECTION("arrayFor")
    {
        std::array<int, 5> arr;
        int count = 0;
        arrayFor(arr, k) { count++; }
        REQUIRE(count == 5);
    }

    SECTION("inFor")
    {
        std::array<std::array<int, 5>, 3> arr;
        int count = 0;
        inFor(arr, k, i) { count++; }
        REQUIRE(count == 15);
    }
}
