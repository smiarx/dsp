#include "dsp/Enveloppe.h"
#include "dsp/Context.h"
#include "dsp/Signal.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>

TEST_CASE("DoubleRamp", "[dsp][doubleramp]")
{
    SECTION("Reference")
    {
        constexpr size_t kN              = 2;
        constexpr dsp::fData<kN> kTarget = {1.f, 2.f};
        constexpr dsp::fData<kN> kUp     = {10.f, 15.f};
        constexpr dsp::fData<kN> kDown   = {40.f, 45.f};

        dsp::DoubleRamp<kN> ramp;
        ramp.set(kTarget, kUp, kDown);

        const dsp::fData<kN> expect[]{
            {0.1f, 0.133333f},   {0.2f, 0.266667f},   {0.3f, 0.4f},
            {0.4f, 0.533333f},   {0.5f, 0.666667f},   {0.6f, 0.8f},
            {0.7f, 0.933333f},   {0.8f, 1.06667f},    {0.9f, 1.2f},
            {1.0f, 1.33333f},    {0.975f, 1.46667f},  {0.95f, 1.6f},
            {0.925f, 1.73333f},  {0.9f, 1.86667f},    {0.875f, 2.0f},
            {0.85f, 1.95556f},   {0.825f, 1.91111f},  {0.8f, 1.86667f},
            {0.775f, 1.82222f},  {0.75f, 1.77778f},   {0.725f, 1.73333f},
            {0.7f, 1.68889f},    {0.675f, 1.64444f},  {0.65f, 1.6f},
            {0.625f, 1.55556f},  {0.6f, 1.51111f},    {0.575f, 1.46667f},
            {0.55f, 1.42222f},   {0.525f, 1.37778f},  {0.5f, 1.33333f},
            {0.475f, 1.28889f},  {0.45f, 1.24444f},   {0.425f, 1.2f},
            {0.4f, 1.15556f},    {0.375f, 1.11111f},  {0.35f, 1.06667f},
            {0.325f, 1.02222f},  {0.3f, 0.977778f},   {0.275f, 0.933333f},
            {0.25f, 0.888889f},  {0.225f, 0.844444f}, {0.2f, 0.8f},
            {0.175f, 0.755556f}, {0.15f, 0.711111f},  {0.125f, 0.666667f},
            {0.1f, 0.622222f},   {0.075f, 0.577778f}, {0.05f, 0.533333f},
            {0.025f, 0.488889f}, {0.0f, 0.444444f},   {0.0f, 0.4f},
            {0.0f, 0.355556f},   {0.0f, 0.311111f},   {0.0f, 0.266667f},
            {0.0f, 0.222222f},   {0.0f, 0.177778f},   {0.0f, 0.133333f},
            {0.0f, 0.0888889f},  {0.0f, 0.0444444f},  {0.f, 0.f},
        };

        size_t n = 0;
        while (ramp.isRunning()) {
            dsp::fSample<kN> x = ramp.process();
            arrayFor(x, i)
            {
                REQUIRE_THAT(x[i],
                             Catch::Matchers::WithinAbs(expect[n][i], 1e-4));
            }
            ++n;
        }
        --n;

        REQUIRE(n == sizeof(expect) / sizeof(expect[0]));
    }
}
