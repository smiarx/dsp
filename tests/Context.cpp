#include "dsp/Context.h"
#include "dsp/MultiVal.h"
#include <array>
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Context Class", "[context]")
{
    SECTION("Construction and Copying")
    {
        dsp::mfloat<> data;
        dsp::Context ctxt(&data);

        // check type
        static_assert(std::is_same_v<decltype(ctxt)::Type, dsp::mfloat<>>);

        REQUIRE((ctxt.getBlockSize() == 1));
        REQUIRE((ctxt.getData() == &data));

        // Copy construction
        dsp::Context ctxtCopy(ctxt);
        REQUIRE((ctxtCopy.getBlockSize() == 1));
        REQUIRE((ctxt.getData() == &data));
    }

    SECTION("Next Method")
    {
        constexpr auto kN = 8;
        dsp::mfloat<> data[kN];
        dsp::Context ctxt(data, kN);

        // check type
        static_assert(std::is_same_v<decltype(ctxt)::Type, dsp::mfloat<>>);

        REQUIRE((ctxt.getBlockSize() == kN));
        ctxt.next();
        REQUIRE((ctxt.getData() == &data[1]));
    }

    SECTION("Data Method")
    {
        constexpr auto kN = 8;
        dsp::mfloat<> data1[kN];
        dsp::mfloat<> data2[kN];
        dsp::Context ctxt(data1, kN);

        REQUIRE((ctxt.getData() == data1));
        ctxt.setData(data2);
        REQUIRE((ctxt.getData() == data2));
    }

    SECTION("Get Input/Output Method")
    {
        // singleton
        float xf      = 0.435f;
        float testVal = 2.34f;
        dsp::Context ctxtf(&xf);
        REQUIRE((xf == ctxtf.getInput()));
        ctxtf.setOutput(testVal);
        REQUIRE((testVal == *ctxtf.getData()));

        // simd
        dsp::mfloat<4> xmf{1, 2, 3, 4};
        dsp::Context ctxt(&xmf);
        auto in = ctxt.getInput();

        // check type
        static_assert(std::is_same_v<decltype(in), dsp::simd<float, 4>>);

        REQUIRE((in[0] == xmf[0]));
        REQUIRE((in[1] == xmf[1]));
        REQUIRE((in[2] == xmf[2]));
        REQUIRE((in[3] == xmf[3]));

        in = sin(in) * 2.f;
        ctxt.setOutput(in);

        auto *data = ctxt.getData();
        REQUIRE((data[0][0] == in[0]));
        REQUIRE((data[0][1] == in[1]));
        REQUIRE((data[0][2] == in[2]));
        REQUIRE((data[0][3] == in[3]));
    }
}

TEST_CASE("Context Run", "[context-run]")
{
    constexpr auto kN            = 21;
    constexpr auto kK            = 2;
    constexpr auto kMaxBatchSize = DSP_MAX_VEC_SIZE / 4;
    std::array<dsp::mfloat<kK>, kN> data{{
        {1.f, 2.f},   {3.f, 4.f},   {5.f, 6.f},   {7.f, 8.f},   {9.f, 10.f},
        {11.f, 12.f}, {13.f, 14.f}, {15.f, 16.f}, {17.f, 18.f}, {19.f, 20.f},
        {21.f, 22.f}, {23.f, 24.f}, {25.f, 26.f}, {27.f, 28.f}, {29.f, 30.f},
        {31.f, 32.f}, {33.f, 34.f}, {35.f, 36.f}, {37.f, 38.f}, {39.f, 40.f},
        {41.f, 42.f},
    }};
    auto orig = data;
    dsp::Context ctxt(data.data(), kN);
    int count = 0;

    SECTION("Scalar")
    {

        dsp::ContextRun(ctxt, [&count](auto c) {
            auto x = c.getInput();
            if (count == 4) {
                REQUIRE((x[0] == 9.f));
                REQUIRE((x[1] == 10.f));
            }
            x *= 2;
            c.setOutput(x);
            count++;
        });

        REQUIRE((count == kN));
    }

    SECTION("Scalar Macro")
    {
        CTXTRUN(ctxt)
        {
            auto x = ctxt.getInput();
            if (count == 6) {
                REQUIRE((x[0] == 13.f));
                REQUIRE((x[1] == 14.f));
            }
            x *= 2;
            ctxt.setOutput(x);
            count++;
        };

        REQUIRE((count == kN));
    }

    SECTION("Vectorized")
    {
        dsp::ContextRun(ctxt.vec(), [&](auto c) {
            constexpr auto kCountCheck = 2;
            constexpr auto kUseVec     = decltype(c)::kUseVec;
            constexpr auto kBatchSize  = kUseVec ? kMaxBatchSize : kK;
            auto x                     = c.getInput();

            if constexpr (kUseVec)
                static_assert(
                    std::is_same_v<decltype(x), dsp::simd<float, kBatchSize>>);

            if (count == kCountCheck) {
                for (int i = 0; i < kBatchSize; ++i) {
                    REQUIRE(
                        (x[i] ==
                         data[(kCountCheck * kBatchSize + i) / kK][i % kK]));
                }
            }

            x *= 2;
            c.setOutput(x);
            count++;
        });

        constexpr auto kNumElements = kN * kK;
        REQUIRE((count == kNumElements / kMaxBatchSize +
                              (kNumElements % kMaxBatchSize) / kK));
    }

    SECTION("Vectorized Macro")
    {
        CTXTRUNVEC(ctxt)
        {
            constexpr auto kCountCheck = 1;
            constexpr auto kUseVec     = decltype(ctxt)::kUseVec;
            constexpr auto kBatchSize  = kUseVec ? kMaxBatchSize : kK;
            auto x                     = ctxt.getInput();

            if constexpr (kUseVec)
                static_assert(
                    std::is_same_v<decltype(x), dsp::simd<float, kBatchSize>>);

            if (count == kCountCheck) {
                for (int i = 0; i < kBatchSize; ++i) {
                    REQUIRE(
                        (x[i] ==
                         data[(kCountCheck * kBatchSize + i) / kK][i % kK]));
                }
            }

            x *= 2;
            ctxt.setOutput(x);
            count++;
        };

        constexpr auto kNumElements = kN * kK;
        REQUIRE((count == kNumElements / kMaxBatchSize +
                              (kNumElements % kMaxBatchSize) / kK));
    }

    for (int i = 0; i < kN; ++i) {
        for (int j = 0; j < kK; ++j) {
            REQUIRE((data[i][j] == orig[i][j] * 2));
        }
    }
}
