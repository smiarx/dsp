#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include "processors/tapedelay/TapeDelay.h"

TEST_CASE("TapeDelay")
{
    SECTION("Benchmark")
    {
        constexpr auto sampleRate = 48000.f;
        constexpr auto blockSize  = 512;

        float in[processors::TapeDelay::N][blockSize];
        float *ins[processors::TapeDelay::N] = {in[0], in[1]};

        processors::TapeDelay tapedelay;
        tapedelay.prepare(sampleRate, blockSize);
        tapedelay.update(0.2, 0.8, 3000, 10, -40, 0.3,
                         processors::TapeDelay::Mode::Normal, 0.4, blockSize);

        BENCHMARK("TapeDelay") { tapedelay.process(ins, ins, blockSize); };

        tapedelay.free();
    }
}
