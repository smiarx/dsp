#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <dsp/core/Context.hpp>
#include <dsp/filters/Notch.hpp>
#include <dsp/generators/Noise.hpp>

//================ Adaptive Notch Filter =========================//

TEST_CASE("RLSANF", "[dsp][rlsanf]")
{
    using T = double;
    dsp::AdaptiveNotchFilter<T> pitchtracker{0.98, 0.995};
    decltype(pitchtracker)::State state{};

    T f0 = GENERATE(0.01, 0.02, 0.1, 0.3);
    T x;

    SECTION("Pure sine")
    {
        for (size_t n = 0; n < 1000; ++n) {
            x = std::sin(dsp::constants<T>::pi * T(n) * f0);
            pitchtracker.process(dsp::Context(&x), state);
        }
        // after convergence we should have the right frequency
        REQUIRE_THAT(state.getFreq(), Catch::Matchers::WithinAbs(f0, 1e-5));
    }
    SECTION("Double sine")
    {
        for (size_t n = 0; n < 1000; ++n) {
            x = std::sin(dsp::constants<T>::pi * T(n) * f0) +
                0.3 * std::sin(dsp::constants<T>::pi * T(n) * 2 * f0);
            pitchtracker.process(dsp::Context(&x), state);
        }
        // after convergence we should have the right frequency
        REQUIRE_THAT(state.getFreq(), Catch::Matchers::WithinAbs(f0, 1e-3));
    }
    SECTION("Sine & noise")
    {
        dsp::Noise<T> noise;
        for (size_t n = 0; n < 1000; ++n) {
            x = std::sin(dsp::constants<T>::pi * T(n) * f0) +
                noise.process() * 0.1;
            pitchtracker.process(dsp::Context(&x), state);
        }
        // after convergence we should have the right frequency
        REQUIRE_THAT(state.getFreq(), Catch::Matchers::WithinAbs(f0, 1e-4));
    }
}
