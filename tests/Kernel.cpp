#include "dsp/Kernel.h"
#include "dsp/Buffer.h"
#include "dsp/Delay.h"
#include "dsp/Windows.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch::Matchers;

template <typename T> static void testIR()
{
    constexpr auto kN = 25;

    constexpr auto kKernelSize = 4;
    using ft                   = T;
    using Kernel = dsp::kernels::Sinc<kKernelSize, dsp::windows::Kaiser<140>>;
    using Tap    = dsp::TapKernel<ft, Kernel, 256>;

    ft data[kN];
    dsp::DelayLine<121> delayline;

    dsp::Buffer<ft, nextTo(delayline)> buffer;
    ft bufferData[decltype(buffer)::kSize]{};
    buffer.setData(bufferData);

    dsp::BufferContext ctxt(data, kN, buffer);

    Tap tap;

    int n                   = 0;
    dsp::baseType<ft> delay = 8.34523;
    int idelay              = static_cast<int>(delay);
    auto fdelay             = delay - idelay;

    tap.setDelay(delay);
    CTXTRUN(ctxt)
    {
        auto x = tap.read(ctxt, delayline);

        ft expect = 0;
        ft sum    = 0;
        for (int i = -kKernelSize; i < kKernelSize; ++i) {
            auto kernel = Kernel::generate(-i - fdelay);
            sum += kernel;
            if (n - idelay + i == 0) {
                expect = dsp::load(expect) + kernel;
            }
        }
        expect /= sum;

        for (size_t i = 0; i < dsp::kTypeWidth<ft>; ++i)
            REQUIRE_THAT(dsp::get(x, i), WithinAbs(dsp::get(expect, i), 1e-5f));
        delayline.write(ctxt, n == 0);
        ++n;
    };
}

template <typename T> static void testDC()
{
    constexpr auto kN = 25;

    constexpr auto kKernelSize = 4;
    using ft                   = T;
    using Kernel = dsp::kernels::Sinc<kKernelSize, dsp::windows::Kaiser<140>>;
    using Tap    = dsp::TapKernel<ft, Kernel, 256>;

    ft data[kN];
    dsp::DelayLine<121> delayline;

    dsp::Buffer<ft, nextTo(delayline)> buffer;
    ft bufferData[decltype(buffer)::kSize]{};
    buffer.setData(bufferData);

    dsp::BufferContext ctxt(data, kN, buffer);

    Tap tap;

    dsp::baseType<ft> delay = 10.38523;
    ft dc                   = 0.5;

    tap.setDelay(delay);
    int n = 0;
    CTXTRUN(ctxt)
    {
        if (n > 2 * kKernelSize + delay) {
            auto x = tap.read(ctxt, delayline);
            for (size_t i = 0; i < dsp::kTypeWidth<ft>; ++i)
                REQUIRE_THAT(dsp::get(x, i), WithinAbs(dsp::get(dc, i), 1e-5f));
        }
        delayline.write(ctxt, dc);
        ++n;
    };
}

template <typename T> static void testBandLimited()
{
    constexpr auto kN = 25;

    constexpr auto kKernelSize = 4;
    using ft                   = T;
    using Kernel = dsp::kernels::Sinc<kKernelSize, dsp::windows::Kaiser<140>>;
    using Tap    = dsp::TapKernel<ft, Kernel, 256>;

    ft data[kN];
    dsp::DelayLine<121> delayline;

    dsp::Buffer<ft, nextTo(delayline)> buffer;
    ft bufferData[decltype(buffer)::kSize]{};
    buffer.setData(bufferData);

    dsp::BufferContext ctxt(data, kN, buffer);

    Tap tap;

    dsp::baseType<ft> delay = 10.38523;

    tap.setDelay(delay);
    int n       = 0;
    float scale = 2;
    SECTION("Nyquist gain is 0")
    {
        CTXTRUN(ctxt)
        {
            auto x = tap.read(ctxt, delayline, scale);
            for (size_t i = 0; i < dsp::kTypeWidth<ft>; ++i)
                REQUIRE_THAT(dsp::get(x, i), WithinAbs(0.0, 1e-5));

            delayline.write(ctxt, dsp::sin(dsp::constants<T>::pi * n));
            ++n;
        };
    }

    SECTION("Reconstruction")
    {
        ft f0 = 0.0191;
        CTXTRUN(ctxt)
        {
            if (n > 2 * kKernelSize + delay) {
                auto x = tap.read(ctxt, delayline, scale);
                auto e = dsp::sin(dsp::constants<T>::pi * (n - delay) *
                                  dsp::load(f0));
                for (size_t i = 0; i < dsp::kTypeWidth<ft>; ++i)
                    REQUIRE_THAT(dsp::get(x, i),
                                 WithinAbs(dsp::get(e, i), 1e-4));
            }

            delayline.write(
                ctxt, dsp::sin(dsp::constants<T>::pi * n * dsp::load(f0)));
            ++n;
        };
    }
}

TEST_CASE("Kernel")
{
    SECTION("impulse response")
    {
        testIR<float>();
        testIR<dsp::mfloat<2>>();
        testIR<dsp::mdouble<2>>();
        testIR<dsp::mfloat<4>>();
    }

    SECTION("DC")
    {
        testDC<float>();
        testDC<dsp::mfloat<2>>();
        testDC<dsp::mdouble<2>>();
        testDC<dsp::mfloat<4>>();
    }

    SECTION("BandLimited")
    {
        testBandLimited<float>();
        testBandLimited<dsp::mfloat<2>>();
        testBandLimited<dsp::mdouble<2>>();
        testBandLimited<dsp::mfloat<4>>();
    }
}
