#include "dsp/Delay.h"
#include "dsp/Buffer.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <iostream>

using namespace Catch::Matchers;

template <typename V>
static inline void checkValue(V delay, int i, float value, float offset = 0.f)
{
    if (i < (int)delay) REQUIRE_THAT(value, WithinRel(0.));
    else if (i > (int)delay)
        REQUIRE_THAT(value,
                     WithinRel(float(i) - static_cast<float>(delay) + offset));
}
template <typename V>
static inline void checkValue(V delay, int i, double value, double offset = 0.)
{
    if (i < (int)delay) REQUIRE_THAT(value, WithinRel(0.));
    else if (i > (int)delay)
        REQUIRE_THAT(
            value, WithinRel(double(i) - static_cast<double>(delay) + offset));
}
template <typename V, typename T>
static inline void checkValue(V delay, int &i, T value, float offset = 0.f)
{
    for (size_t j = 0; j < decltype(value)::kWidth; ++j) {
        checkValue(delay, i + static_cast<int>(j), value[j], offset);
    }
}

TEST_CASE("DelayLine", "[dsp][delayline]")
{
    using ft = float;

    constexpr auto kDelay1    = 10;
    constexpr auto kDelay2    = 19;
    constexpr auto kN         = 50;
    constexpr auto kBlockSize = 13;

    dsp::DelayLine<kDelay1 + 5> delay1;
    dsp::DelayLine<kDelay2 + 7, nextTo(delay1)> delay2;
    dsp::Buffer<ft, nextTo(delay2)> buffer;

    ft bufdata[decltype(buffer)::kSize]{};
    buffer.setData(bufdata);

    ft data[kN];
    for (size_t i = 0; i < kN; ++i) {
        data[i] = static_cast<ft>(i);
    }

    SECTION("Multiple Delays")
    {
        int i = 0;

        for (size_t j = 0; j < kN / kBlockSize; ++j) {
            dsp::BufferContext ctxt(data + kBlockSize * j, kBlockSize, buffer);

            CTXTRUN(ctxt)
            {
                auto x = ctxt.getInput();

                auto v1 = delay1.read(ctxt, kDelay1);
                auto v2 = delay2.read(ctxt, kDelay2);

                delay1.write(ctxt, x);
                delay2.write(ctxt, x + 10.f);

                checkValue(kDelay1, i, v1);
                checkValue(kDelay2, i, v2, 10.f);
                ++i;
            };
            buffer.nextBlock(ctxt);
        }
    }

    SECTION("Multiple Delays Vectorized")
    {
        int i = 0;

        for (size_t j = 0; j < kN / kBlockSize; ++j) {
            dsp::BufferContext ctxt(data + kBlockSize * j, kBlockSize, buffer);

            CTXTRUNVEC(ctxt)
            {
                auto x  = ctxt.getInput();
                auto v1 = delay1.read(ctxt, kDelay1);
                auto v2 = delay2.read(ctxt, kDelay2);

                delay1.write(ctxt, x);
                delay2.write(ctxt, x + 10.f);

                checkValue(kDelay1, i, v1);
                checkValue(kDelay2, i, v2, 10.f);
                if constexpr (std::is_same_v<ft, decltype(v1)>) ++i;
                else
                    i += int(decltype(v1)::kWidth);
            };
            buffer.nextBlock(ctxt);
        }
    }
}

TEST_CASE("CopyDelayLine", "[dsp][copydelayline]")
{

    using namespace Catch::Matchers;

    using ft = float;

    constexpr auto kDelay     = 13;
    constexpr auto kBlockSize = 23;
    constexpr auto kN         = kBlockSize * 3;

    dsp::CopyDelayLine<ft, kDelay + 10> delay;
    ft data[kBlockSize];

    SECTION("Scalar")
    {

        int i = 0;
        for (size_t j = 0; j < kN / kBlockSize; ++j) {

            for (size_t k = 0; k < kBlockSize; ++k) {
                data[k] = static_cast<ft>(k + j * kBlockSize);
            }

            dsp::Context ctxt(data, kBlockSize);

            CTXTRUN(ctxt)
            {
                auto x = ctxt.getInput();

                auto v = delay.read(ctxt, kDelay);
                delay.write(ctxt, x);

                checkValue(kDelay, i, v);
                ++i;
            };
        }
    }

    SECTION("Vector")
    {

        int i = 0;
        for (size_t j = 0; j < kN / kBlockSize; ++j) {

            for (size_t k = 0; k < kBlockSize; ++k) {
                data[k] = static_cast<ft>(k + j * kBlockSize);
            }

            dsp::Context ctxt(data, kBlockSize);

            CTXTRUNVEC(ctxt)
            {
                auto x = ctxt.getInput();

                auto v = delay.read(ctxt, kDelay);
                delay.write(ctxt, x);

                checkValue(kDelay, i, v);
                if constexpr (std::is_same_v<ft, decltype(x)>) ++i;
                else
                    i += static_cast<int>(decltype(x)::kWidth);
            };
        }
    }
}

TEST_CASE("Tap", "[dsp][delay][tap]")
{

    SECTION("Scalar")
    {
        using namespace Catch::Matchers;

        using ft = float;

        static constexpr auto kDelayMax  = 13;
        static constexpr auto kBlockSize = 23;
        static constexpr auto kN         = kBlockSize * 3;

        ft data[kBlockSize];

        dsp::DelayLine<kDelayMax> delay;
        dsp::Buffer<ft, nextTo(delay)> buffer;

        ft bufdata[decltype(buffer)::kSize]{};
        buffer.setData(bufdata);

        int i = 0;

        for (size_t j = 0; j < kN / kBlockSize; ++j) {

            for (size_t k = 0; k < kBlockSize; ++k) {
                data[k] = static_cast<ft>(k + j * kBlockSize);
            }

            dsp::BufferContext ctxt(data, kBlockSize, buffer);

            CTXTRUN(ctxt)
            {
                auto x = ctxt.getInput();

                // TapTail
                {
                    auto v = dsp::TapTail{}.read(ctxt, delay);
                    checkValue(kDelayMax, i, v);
                }
                // TapFix
                {
                    constexpr auto kD = 5;
                    auto v            = dsp::TapFix<kD>{}.read(ctxt, delay);
                    checkValue(kD, i, v);
                }
                // TapNoInterp
                {
                    auto d = GENERATE(take(1, random(1, kDelayMax)));
                    auto v = dsp::TapNoInterp<ft>(d).read(ctxt, delay);
                    checkValue(d, i, v);
                }
                // TapLin
                {
                    auto d = GENERATE(take(1, random(1.f, ft(kDelayMax - 1))));
                    auto v = dsp::TapLin<ft>(d).read(ctxt, delay);
                    checkValue(d, i, v);
                }
                // TapCubic
                {
                    auto d = GENERATE(take(1, random(2.f, ft(kDelayMax - 2))));
                    auto v = dsp::TapCubic<ft>(d).read(ctxt, delay);

                    if (ft(i) < d - 1.f) {
                        REQUIRE_THAT(v, WithinRel(ft(0)));
                    } else {
                        ft x0     = ft(i) - ft(static_cast<int>(d));
                        ft x1     = std::max(ft(0), x0 - 1);
                        ft x2     = std::max(ft(0), x1 - 1);
                        ft xm1    = x0 + 1;
                        ft fd     = d - ft(static_cast<int>(d));
                        ft expect = dsp::hermite(xm1, x0, x1, x2, fd);

                        REQUIRE_THAT(v, WithinRel(expect));
                    }
                }

                delay.write(ctxt, x);
                ++i;
            };
            buffer.nextBlock(ctxt);
        }
    }

    SECTION("Vector")
    {
        using namespace Catch::Matchers;

        using ft = float;

        static constexpr auto kDelayMax  = 23;
        static constexpr auto kBlockSize = 23;
        static constexpr auto kN         = kBlockSize * 3;

        ft data[kBlockSize];

        dsp::DelayLine<kDelayMax> delay;
        dsp::Buffer<ft, nextTo(delay)> buffer;

        ft bufdata[decltype(buffer)::kSize]{};
        buffer.setData(bufdata);

        int i = 0;

        for (size_t j = 0; j < kN / kBlockSize; ++j) {

            for (size_t k = 0; k < kBlockSize; ++k) {
                data[k] = static_cast<ft>(k + j * kBlockSize);
            }

            dsp::BufferContext ctxt(data, kBlockSize, buffer);

            CTXTRUNVEC(ctxt)
            {
                constexpr auto kMinDelay = DSP_MAX_VEC_SIZE / 4;
                auto x                   = ctxt.getInput();

                // TapTail
                {
                    auto v = dsp::TapTail{}.read(ctxt, delay);
                    checkValue(kDelayMax, i, v);
                }
                // TapFix
                {
                    constexpr auto kD = 10;
                    auto v            = dsp::TapFix<kD>{}.read(ctxt, delay);
                    checkValue(kD, i, v);
                }
                // TapNoInterp
                {
                    auto d = GENERATE(
                        take(1, random(kMinDelay, kDelayMax - kMinDelay)));
                    auto v = dsp::TapNoInterp<ft>(d).read(ctxt, delay);
                    static_assert(
                        std::is_same_v<decltype(v), float> ||
                        std::is_same_v<decltype(v), dsp::mfloat<>::simdtype>);
                    checkValue(d, i, v);
                }
                // TapLin
                {
                    auto d = GENERATE(take(
                        1, random(ft(kMinDelay), ft(kDelayMax - kMinDelay))));
                    auto v = dsp::TapLin<ft>(d).read(ctxt, delay);
                    checkValue(d, i, v);
                }

                delay.write(ctxt, x);
                if constexpr (std::is_same_v<ft, decltype(x)>) ++i;
                else
                    i += static_cast<int>(decltype(x)::kWidth);
            };
            buffer.nextBlock(ctxt);
        }
    }

    SECTION("Multi Channel")
    {
        using namespace Catch::Matchers;

#ifdef DSP_SIMD_DOUBLE
        using ft    = dsp::mdouble<2>;
        using fscal = double;
#else
        using ft    = dsp::mfloat<2>;
        using fscal = float;
#endif

        constexpr auto kDelayMax  = 13;
        constexpr auto kBlockSize = 29;
        constexpr auto kNBlocks   = 2;

        ft data[kBlockSize];

        dsp::DelayLine<kDelayMax> delay;
        dsp::Buffer<ft, nextTo(delay)> buffer;

        ft bufdata[decltype(buffer)::kSize]{};
        buffer.setData(bufdata);

        int i     = 0;
        fscal off = 6.;

        for (size_t j = 0; j < kNBlocks; ++j) {

            for (size_t k = 0; k < kBlockSize; ++k) {
                data[k][0] = static_cast<fscal>(k + j * kBlockSize);
                data[k][1] = static_cast<fscal>(k + j * kBlockSize) + off;
            }

            dsp::BufferContext ctxt(data, kBlockSize, buffer);

            CTXTRUN(ctxt)
            {
                auto x = ctxt.getInput();

                // TapFix
                {
                    constexpr int kD[2] = {5, 9};
                    auto v = dsp::TapFix<kD[0], kD[1]>{}.read(ctxt, delay);
                    checkValue(kD[0], i, v[0]);
                    checkValue(kD[1], i, v[1], off);
                }
                // TapNoInterp
                {
                    dsp::mint<2> d = {GENERATE(take(1, random(1, kDelayMax))),
                                      GENERATE(take(1, random(1, kDelayMax)))};
                    auto v         = dsp::TapNoInterp<ft>(d).read(ctxt, delay);
                    checkValue(d[0], i, v[0]);
                    checkValue(d[1], i, v[1], off);
                }
                // TapLin
                {
                    ft d   = {GENERATE(take(
                                1, random(fscal(1.),
                                            static_cast<fscal>(kDelayMax - 1)))),
                              GENERATE(take(
                                1, random(fscal(1.),
                                            static_cast<fscal>(kDelayMax - 1))))};
                    auto v = dsp::TapLin<ft>(d).read(ctxt, delay);
                    checkValue(d[0], i, v[0]);
                    checkValue(d[1], i, v[1], off);
                }
                // TapCubic
                {
                    ft d   = {GENERATE(take(
                                1, random(fscal(2.),
                                            static_cast<fscal>(kDelayMax - 2)))),
                              GENERATE(take(
                                1, random(fscal(2.),
                                            static_cast<fscal>(kDelayMax - 2))))};
                    auto v = dsp::TapCubic<ft>(d).read(ctxt, delay);

                    for (size_t j = 0; j < 2; ++j) {
                        fscal x0  = i - static_cast<int>(d[j]) + off;
                        fscal x1  = x0 - 1;
                        fscal x2  = x1 - 1;
                        fscal xm1 = x0 + 1;
                        if (j == 1) {
                            // set offset
                            x0  = x0 < off ? 0. : x0;
                            x1  = x1 < off ? 0. : x1;
                            x2  = x2 < off ? 0. : x2;
                            xm1 = xm1 < off ? 0. : xm1;
                        }
                        auto fd     = d[j] - static_cast<int>(d[j]);
                        auto expect = dsp::hermite(xm1, x0, x1, x2, fd);
                        if (j == 1) REQUIRE_THAT(v[j], WithinRel(expect));
                    }
                }

                delay.write(ctxt, x);
                ++i;
            };
            buffer.nextBlock(ctxt);
        }
    }
}
