#include "dsp/RMS.h"
#include "dsp/Stack.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

TEST_CASE("RMS", "[dsp][rms]")
{
    using namespace Catch::Generators;
    using namespace Catch::Matchers;

    SECTION("DC")
    {
        constexpr size_t N       = 1;
        constexpr size_t RMSSize = 256;
        constexpr size_t Size    = 32;
        constexpr size_t NBlock  = RMSSize / Size;
        constexpr size_t NValues = 4;

        dsp::RMS<N, RMSSize> rms;
        dsp::Stack<N, 1> stack;
        dsp::fSample<N> x[Size];
        float values[NValues] = {GENERATE(take(4, random(-10.f, 10.f))),
                                 GENERATE(take(4, random(-10.f, 10.f))),
                                 GENERATE(take(4, random(-10.f, 10.f))),
                                 GENERATE(take(4, random(-10.f, 10.f)))};

        for (size_t t = 0; t < NValues; ++t) {
            for (size_t n = 0; n < Size; ++n) {
                x[n][0] = values[t];
            }

            for (size_t k = 0; k < NBlock; ++k) {
                rms.processBlock(dsp::Context(x, Size), stack);
            }

            REQUIRE_THAT(stack.get()[0], WithinAbs(std::abs(values[t]), 1e-4));
        }
    }

    SECTION("sine")
    {
        constexpr size_t N       = 1;
        constexpr size_t RMSSize = 256;
        constexpr size_t Size    = 32;
        constexpr size_t NBlock  = RMSSize / Size;

        dsp::RMS<N, RMSSize> rms;
        dsp::Stack<N, 1> stack;
        dsp::fSample<N> x[Size];
        float freq = GENERATE(take(4, random(4.f / RMSSize, 1.f)));

        for (size_t t = 0; t < 2; ++t) {
            for (size_t k = 0; k < NBlock; ++k) {
                for (size_t n = 0; n < Size; ++n) {
                    x[n][0] = sin(dsp::constants<float>::pi * freq * n);
                }

                rms.processBlock(dsp::Context(x, Size), stack);
            }

            REQUIRE_THAT(stack.get()[0], WithinAbs(std::abs(M_SQRT1_2f), 1e-1));
        }
    }

    SECTION("Overlap")
    {
        constexpr size_t N        = 1;
        constexpr size_t RMSSize  = 256;
        constexpr size_t Shift    = 16;
        constexpr size_t Overlap  = RMSSize - Shift;
        constexpr size_t Size     = 8;
        constexpr size_t NBlock   = Shift / Size;
        constexpr size_t NOverlap = RMSSize / Shift;
        constexpr size_t NValues  = 4;

        dsp::RMS<N, RMSSize, Overlap> rms;
        dsp::Stack<N, 1> stack;
        dsp::fSample<N> x[Size];

        for (size_t n = 0; n < Size; ++n) {
            x[n][0] = 1.f;
        }

        for (size_t j = 0; j < NOverlap; ++j) {
            for (size_t k = 0; k < NBlock; ++k) {
                rms.processBlock(dsp::Context(x, Size), stack);
            }
            REQUIRE_THAT(stack.get()[0],
                         WithinAbs(sqrtf((1.f + j) * Shift / RMSSize), 1e-6));
        }

        for (size_t n = 0; n < Size; ++n) {
            x[n][0] = 0.f;
        }

        for (size_t j = 0; j < NOverlap; ++j) {
            for (size_t k = 0; k < NBlock; ++k) {
                rms.processBlock(dsp::Context(x, Size), stack);
            }
            REQUIRE_THAT(
                stack.get()[0],
                WithinAbs(sqrtf((NOverlap - 1.f - j) * Shift / RMSSize), 1e-6));
        }
    }
}
