#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include "processors/tapedelay/TapeDelay.h"

TEST_CASE("TapeDelay")
{
    SECTION("Benchmark")
    {
        constexpr auto kSampleRate = 48000.f;
        constexpr auto kBlockSize  = 512;

        float in[processors::TapeDelay::kN][kBlockSize];
        float *ins[processors::TapeDelay::kN] = {in[0], in[1]};

        processors::TapeDelay tapedelay;
        tapedelay.prepare(kSampleRate, kBlockSize);
        tapedelay.update(0.2f, 0.8f, 3000.f, 10.f, -40.f, 0.3f,
                         processors::TapeDelay::Mode::kNormal, 0.4f,
                         kBlockSize);

        BENCHMARK("TapeDelay") { tapedelay.process(ins, ins, kBlockSize); };

        tapedelay.free();
    }
}
