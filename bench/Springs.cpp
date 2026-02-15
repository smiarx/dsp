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
        springs.update(0.5f, 4500.f, 0.05f, 2.3f, 0.3f, 0.4f, 0.4f, 1.f, 0.3f,
                       kBlockSize);
        BENCHMARK("Springs") { springs.process(ins, ins, kBlockSize); };
        springs.free();
    }
}
