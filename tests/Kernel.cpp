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
#include <iostream>

using namespace Catch::Matchers;

template <typename T> static void testKernel()
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
        for (int i = -kKernelSize; i < kKernelSize; ++i) {
            if (n - idelay + i == 0)
                expect = dsp::load(expect) + Kernel::generate(-i - fdelay);
        }
        for (size_t i = 0; i < dsp::kTypeWidth<ft>; ++i)
            REQUIRE_THAT(dsp::get(x, i), WithinAbs(dsp::get(expect, i), 1e-5f));
        delayline.write(ctxt, n == 0);
        ++n;
    };
}

TEST_CASE("Kernel")
{
    testKernel<float>();
    testKernel<dsp::mfloat<2>>();
    testKernel<dsp::mdouble<2>>();
    testKernel<dsp::mfloat<4>>();
}
