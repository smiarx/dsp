#include "dsp/RMS.h"
#include "dsp/Context.h"
#include "dsp/FastMath.h"
#include "dsp/Stack.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstddef>

TEST_CASE("RMS", "[dsp][rms]")
{
    using namespace Catch::Generators;
    using namespace Catch::Matchers;

    SECTION("DC")
    {
        using ft                  = float;
        constexpr size_t kRmsSize = 256;
        constexpr size_t kSize    = 32;
        constexpr size_t kNBlock  = kRmsSize / kSize;
        constexpr size_t kNValues = 4;

        dsp::RMS<ft, kRmsSize> rms;
        dsp::Stack<ft, 1> stack;
        ft x[kSize];
        float values[kNValues] = {GENERATE(take(4, random(-10.f, 10.f))),
                                  GENERATE(take(4, random(-10.f, 10.f))),
                                  GENERATE(take(4, random(-10.f, 10.f))),
                                  GENERATE(take(4, random(-10.f, 10.f)))};

        for (float value : values) {
            for (auto &xn : x) {
                xn = value;
            }

            for (size_t k = 0; k < kNBlock; ++k) {
                rms.processBlock(dsp::Context(x, kSize), stack);
            }

            REQUIRE_THAT(stack.get(), WithinAbs(std::abs(value), 1e-4));
        }
    }

    SECTION("sine")
    {
        using ft                  = float;
        constexpr size_t kRmsSize = 256;
        constexpr size_t kSize    = 32;
        constexpr size_t kNBlock  = kRmsSize / kSize;

        dsp::RMS<ft, kRmsSize> rms;
        dsp::Stack<ft, 1> stack;
        ft x[kSize];
        float freq = GENERATE(take(4, random(4.f / kRmsSize, 1.f)));

        for (size_t t = 0; t < 2; ++t) {
            for (size_t k = 0; k < kNBlock; ++k) {
                for (size_t n = 0; n < kSize; ++n) {
                    x[n] = std::sin(dsp::constants<float>::pi * freq * n);
                }

                rms.processBlock(dsp::Context(x, kSize), stack);
            }

            REQUIRE_THAT(
                stack.get(),
                WithinAbs(std::abs(dsp::constants<float>::sqrt1_2), 1e-1));
        }
    }

    SECTION("Overlap")
    {
        using ft                   = float;
        constexpr size_t kRmsSize  = 256;
        constexpr size_t kShift    = 16;
        constexpr size_t kOverlap  = kRmsSize - kShift;
        constexpr size_t kSize     = 8;
        constexpr size_t kNBlock   = kShift / kSize;
        constexpr size_t kNOverlap = kRmsSize / kShift;

        dsp::RMS<ft, kRmsSize, kOverlap> rms;
        dsp::Stack<ft, 1> stack;
        ft x[kSize];

        for (auto &xn : x) {
            xn = 1.f;
        }

        for (size_t j = 0; j < kNOverlap; ++j) {
            for (size_t k = 0; k < kNBlock; ++k) {
                rms.processBlock(dsp::Context(x, kSize), stack);
            }
            REQUIRE_THAT(stack.get(),
                         WithinAbs(sqrtf((1.f + j) * kShift / kRmsSize), 1e-6));
        }

        for (auto &xn : x) {
            xn = 0.f;
        }

        for (size_t j = 0; j < kNOverlap; ++j) {
            for (size_t k = 0; k < kNBlock; ++k) {
                rms.processBlock(dsp::Context(x, kSize), stack);
            }
            REQUIRE_THAT(stack.get(), WithinAbs(sqrtf((kNOverlap - 1.f - j) *
                                                      kShift / kRmsSize),
                                                1e-6));
        }
    }
}
