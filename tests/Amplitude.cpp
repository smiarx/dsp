#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dsp/core/Context.hpp>
#include <dsp/filters/Amplitude.hpp>

using namespace Catch::Matchers;

TEST_CASE("Amplitude")
{
    SECTION("float")
    {
        dsp::Amplitude<float> amp{};
        decltype(amp)::State state{};

        constexpr auto kNUp    = 30;
        constexpr auto kNDown  = 100;
        constexpr auto kTarget = 0.01f;

        amp.setRate(kTarget, static_cast<float>(kNUp),
                    static_cast<float>(kNDown));

        for (int i = 0; i < kNUp; ++i) {
            auto one = 1.f;
            amp.process(dsp::Context(&one), state);
        }
        REQUIRE_THAT(state, WithinAbs(1.f - kTarget, 1e-6));

        state = 1.f;
        for (int i = 0; i < kNDown; ++i) {
            auto zero = 0.f;
            amp.process(dsp::Context(&zero), state);
        }
        REQUIRE_THAT(state, WithinAbs(kTarget, 1e-6));
    }

    SECTION("float4")
    {
        using type = dsp::mfloat4;
        dsp::Amplitude<type> amp{};
        decltype(amp)::State state{};

        constexpr type kNUp    = {20.f, 30.f, 33.f, 44.f};
        constexpr type kNDown  = {40.f, 80.f, 100.f, 300.f};
        constexpr type kTarget = {0.01f, 0.001f, 0.05f, 0.2f};

        amp.setRate(kTarget, kNUp, kNDown);

        constexpr auto kMaxUp   = static_cast<int>(kNUp[3]);
        constexpr auto kMaxDown = static_cast<int>(kNDown[3]);

        size_t i = 0;
        for (int n = 0; n < kMaxUp; ++n) {
            if (n == static_cast<int>(kNUp[i])) {
                REQUIRE_THAT(state[i], WithinAbs(1.f - kTarget[i], 1e-6));
                ++i;
            }
            type one = 1.f;
            amp.process(dsp::Context(&one), state);
        }
        REQUIRE_THAT(state[i], WithinAbs(1.f - kTarget[i], 1e-6));

        i     = 0;
        state = 1.f;
        for (int n = 0; n < kMaxDown; ++n) {
            if (n == static_cast<int>(kNDown[i])) {
                REQUIRE_THAT(state[i], WithinAbs(kTarget[i], 1e-6));
                ++i;
            }
            type zero = 0.f;
            amp.process(dsp::Context(&zero), state);
        }
        REQUIRE_THAT(state[i], WithinAbs(kTarget[i], 1e-6));
    }
}
