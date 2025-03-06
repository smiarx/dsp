#include "dsp/Signal.h"
#include <catch2/catch_test_macros.hpp>
#include <cstdlib> // for std::aligned_alloc and std::free

TEST_CASE("Data Class", "[data]")
{
    SECTION("Data<float, 4>")
    {
        dsp::Data<float, 4> d({1.0f, 2.0f, 3.0f, 4.0f});
        auto v = d.toSIMD();
        float result[4];
        dsp::SIMDtoArray<float, 4, true>(result, v);
        REQUIRE(result[0] == 1.0f);
        REQUIRE(result[1] == 2.0f);
        REQUIRE(result[2] == 3.0f);
        REQUIRE(result[3] == 4.0f);
    }

    SECTION("Data<int, 4>")
    {
        dsp::Data<int, 4> d({1, 2, 3, 4});
        auto v = d.toSIMD();
        int result[4];
        dsp::SIMDtoArray<int, 4, true>(result, v);
        REQUIRE(result[0] == 1);
        REQUIRE(result[1] == 2);
        REQUIRE(result[2] == 3);
        REQUIRE(result[3] == 4);
    }
}

TEST_CASE("Sample Class", "[sample]")
{
    SECTION("Sample<float, 4>")
    {
        dsp::Sample<float, 4> s({1.0f, 2.0f, 3.0f, 4.0f});

        // Convert to scalar signal
        auto &scalar_signal = s.toScalar();
        REQUIRE(scalar_signal[0][0] == 1.0f);
        REQUIRE(scalar_signal[0][1] == 2.0f);
        REQUIRE(scalar_signal[0][2] == 3.0f);
        REQUIRE(scalar_signal[0][3] == 4.0f);

        // Convert to vectorized signal
        auto &vector_signal = s.toVector();
        REQUIRE(vector_signal[0][0] == 1.0f);
        REQUIRE(vector_signal[0][1] == 2.0f);
        REQUIRE(vector_signal[0][2] == 3.0f);
        REQUIRE(vector_signal[0][3] == 4.0f);
    }

    SECTION("Sample<int, 4>")
    {
        dsp::Sample<int, 4> s({1, 2, 3, 4});

        // Convert to scalar signal
        auto &scalar_signal = s.toScalar();
        REQUIRE(scalar_signal[0][0] == 1);
        REQUIRE(scalar_signal[0][1] == 2);
        REQUIRE(scalar_signal[0][2] == 3);
        REQUIRE(scalar_signal[0][3] == 4);

        // Convert to vectorized signal
        auto &vector_signal = s.toVector();
        REQUIRE(vector_signal[0][0] == 1);
        REQUIRE(vector_signal[0][1] == 2);
        REQUIRE(vector_signal[0][2] == 3);
        REQUIRE(vector_signal[0][3] == 4);
    }
}

TEST_CASE("Signal Class", "[signal]")
{
    SECTION("Signal<float, 2, 2>")
    {
        float data[] = {1.0f, 2.0f, 3.0f, 4.0f};
        dsp::Signal<float, 2, 2> s;
        std::copy(data, data + 4, s.data()->data());

        auto v = s.toSIMD();
        float result[4];
        dsp::SIMDtoArray<float, 4, true>(result, v);
        REQUIRE(result[0] == 1.0f);
        REQUIRE(result[1] == 2.0f);
        REQUIRE(result[2] == 3.0f);
        REQUIRE(result[3] == 4.0f);

        // Aligned load and store
        float *aligned_data = (float *)_mm_malloc(4 * sizeof(float), 16);
        std::copy(data, data + 4, aligned_data);
        dsp::Signal<float, 2, 2> s_aligned;
        std::copy(aligned_data, aligned_data + 4, s_aligned.data()->data());
        v = s_aligned.toSIMD();
        dsp::SIMDtoArray<float, 4, true>(aligned_data, v);
        REQUIRE(aligned_data[0] == 1.0f);
        REQUIRE(aligned_data[1] == 2.0f);
        REQUIRE(aligned_data[2] == 3.0f);
        REQUIRE(aligned_data[3] == 4.0f);
        _mm_free(aligned_data);
    }

#if defined(__AVX__)
    SECTION("Signal<float, 4, 2>")
    {
        float data[] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};
        dsp::Signal<float, 4, 2> s;
        std::copy(data, data + 8, s.data()->data());

        auto v = s.toSIMD();
        float result[8];
        SIMDtoArray<float, 8, true>(result, v);
        REQUIRE(result[0] == 1.0f);
        REQUIRE(result[1] == 2.0f);
        REQUIRE(result[2] == 3.0f);
        REQUIRE(result[3] == 4.0f);
        REQUIRE(result[4] == 5.0f);
        REQUIRE(result[5] == 6.0f);
        REQUIRE(result[6] == 7.0f);
        REQUIRE(result[7] == 8.0f);

        // Aligned load and store
        float *aligned_data = (float *)_mm_malloc(8 * sizeof(float), 32);
        std::copy(data, data + 8, aligned_data);
        dsp::Signal<float, 4, 2> s_aligned;
        std::copy(aligned_data, aligned_data + 8, s_aligned.data()->data());
        v = s_aligned.toSIMD();
        SIMDtoArray<float, 8, true>(aligned_data, v);
        REQUIRE(aligned_data[0] == 1.0f);
        REQUIRE(aligned_data[1] == 2.0f);
        REQUIRE(aligned_data[2] == 3.0f);
        REQUIRE(aligned_data[3] == 4.0f);
        REQUIRE(aligned_data[4] == 5.0f);
        REQUIRE(aligned_data[5] == 6.0f);
        REQUIRE(aligned_data[6] == 7.0f);
        REQUIRE(aligned_data[7] == 8.0f);
        _mm_free(aligned_data);
    }

    SECTION("Signal<int, 4, 2>")
    {
        int data[] = {1, 2, 3, 4, 5, 6, 7, 8};
        dsp::Signal<int, 4, 2> s;
        std::copy(data, data + 8, s.data()->data());

        auto v = s.toSIMD();
        int result[8];
        SIMDtoArray<int, 8, true>(result, v);
        REQUIRE(result[0] == 1);
        REQUIRE(result[1] == 2);
        REQUIRE(result[2] == 3);
        REQUIRE(result[3] == 4);
        REQUIRE(result[4] == 5);
        REQUIRE(result[5] == 6);
        REQUIRE(result[6] == 7);
        REQUIRE(result[7] == 8);

        // Aligned load and store
        int *aligned_data = (int *)_mm_malloc(8 * sizeof(int), 32);
        std::copy(data, data + 8, aligned_data);
        dsp::Signal<int, 4, 2> s_aligned;
        std::copy(aligned_data, aligned_data + 8, s_aligned.data()->data());
        v = s_aligned.toSIMD();
        SIMDtoArray<int, 8, true>(aligned_data, v);
        REQUIRE(aligned_data[0] == 1);
        REQUIRE(aligned_data[1] == 2);
        REQUIRE(aligned_data[2] == 3);
        REQUIRE(aligned_data[3] == 4);
        REQUIRE(aligned_data[4] == 5);
        REQUIRE(aligned_data[5] == 6);
        REQUIRE(aligned_data[6] == 7);
        REQUIRE(aligned_data[7] == 8);
        _mm_free(aligned_data);
    }
#endif
}
