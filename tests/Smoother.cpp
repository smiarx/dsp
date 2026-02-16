#include "dsp/Smoother.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>

template <typename T, bool Vec = false, bool VecCtxt = Vec>
static void testControlSmoother()
{
    using bT              = dsp::baseType<T>;
    const auto N          = 12;
    const auto rate       = bT(1) / bT(N);
    const T init          = T(1.);
    const T value         = T(2.);
    const auto eps        = bT(1e-6);
    constexpr auto kWidth = dsp::kTypeWidth<T>;

    dsp::ControlSmoother<T, Vec> smoother{init};

    REQUIRE(!smoother.isActive());
    smoother.set(value, rate);
    REQUIRE(smoother.isActive());

    dsp::Context<T, VecCtxt> ctxt(nullptr, N);

    auto x = smoother.get();
    for (size_t i = 0; i < dsp::kTypeWidth<T>; ++i)
        REQUIRE_THAT(dsp::get(x, i),
                     Catch::Matchers::WithinRel(dsp::get(init, i), eps));

    size_t k = 0;
    CTXTRUN(ctxt)
    {
        auto x = smoother.step(ctxt);

        constexpr auto kIncrSize = decltype(ctxt)::kIncrSize;
        for (size_t j = 0; j < kIncrSize; ++j) {
            auto expect = init + (value - dsp::load(init)) * rate *
                                     static_cast<T>(k + j + 1);

            for (size_t i = 0; i < kWidth; ++i)
                REQUIRE_THAT(
                    dsp::get(x, j * kWidth + i),
                    Catch::Matchers::WithinRel(dsp::get(expect, i), eps));
        }

        k += kIncrSize;
    };

    // last value is target
    x = smoother.get();
    for (size_t i = 0; i < kWidth; ++i)
        REQUIRE_THAT(dsp::get(x, i),
                     Catch::Matchers::WithinRel(dsp::get(value, i), eps));

    // inactive after reset
    smoother.reset();
    REQUIRE(!smoother.isActive());

    // dont do anything if inactive
    smoother.step(ctxt);
    x = smoother.get();
    for (size_t i = 0; i < dsp::kTypeWidth<T>; ++i)
        REQUIRE_THAT(dsp::get(x, i),
                     Catch::Matchers::WithinRel(dsp::get(value, i), eps));

    // dont activate if target doesn't change
    smoother.set(value, rate);
    REQUIRE(!smoother.isActive());
}

TEST_CASE("Control Smoother", "[dsp][smoother][controls]")
{
    SECTION("float") { testControlSmoother<float>(); }
    SECTION("double") { testControlSmoother<double>(); }
    SECTION("float x 2") { testControlSmoother<dsp::mfloat<2>>(); }
    SECTION("double x 2") { testControlSmoother<dsp::mdouble<2>>(); }
    SECTION("float x 4") { testControlSmoother<dsp::mfloat<4>>(); }

    SECTION("float vec") { testControlSmoother<float, true>(); }
    SECTION("float x 2 vec") { testControlSmoother<dsp::mfloat<2>, true>(); }
    SECTION("float x 4 vec") { testControlSmoother<dsp::mfloat<4>, true>(); }

    SECTION("float vec with scalar context")
    {
        testControlSmoother<float, true, false>();
    }
    SECTION("float x 2 vec with scalar context")
    {
        testControlSmoother<dsp::mfloat<2>, true, false>();
    }
    SECTION("float x 4 vec with scalar context")
    {
        testControlSmoother<dsp::mfloat<4>, true, false>();
    }
}
