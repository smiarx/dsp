#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include "../processors/springs/Springs.h"

TEST_CASE("Springs")
{
    SECTION("Benchmark")
    {
        constexpr auto kSampleRate = 48000.f;
        constexpr auto kBlockSize  = 512;

        float in[2][kBlockSize];
        float *ins[2] = {in[0], in[1]};

        processors::Springs springs;
        springs.prepare(kSampleRate, kBlockSize);
        springs.update(0.5, 4500, 0.05, 2.3, 0.3, 0.4, 0.4, 1., 0.3,
                       kBlockSize);
        BENCHMARK("Springs") { springs.process(ins, ins, kBlockSize); };
        springs.free();
    }
}
