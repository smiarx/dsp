#include "dsp/SIMD.h"
#include <catch2/catch_test_macros.hpp>
#include <cstdlib> // for std::aligned_alloc and std::free

TEST_CASE("SIMD Loading and Storing", "[simd]")
{
    SECTION("SIMD<int, 1>")
    {
        int data[] = {1};
        auto v     = dsp::arrayToSIMD<int, 1>(data);
        REQUIRE(v == 1);

        int result[1];
        dsp::SIMDtoArray<int, 1>(result, v);
        REQUIRE(result[0] == 1);
    }

    SECTION("SIMD<int, 2> SSE")
    {
#if defined(__SSE__)
        int data[] = {1, 2};
        auto v     = dsp::arrayToSIMD<int, 2>(data);
        int result[2];
        dsp::SIMDtoArray<int, 2>(result, v);
        REQUIRE(result[0] == 1);
        REQUIRE(result[1] == 2);

        // Aligned load and store
        int *alignedData = (int *)_mm_malloc(2 * sizeof(int), 16);
        alignedData[0]   = 1;
        alignedData[1]   = 2;
        v                = dsp::arrayToSIMD<int, 2, true>(alignedData);
        dsp::SIMDtoArray<int, 2, true>(alignedData, v);
        REQUIRE(alignedData[0] == 1);
        REQUIRE(alignedData[1] == 2);
        _mm_free(alignedData);
#endif
    }

    SECTION("SIMD<int, 4> SSE")
    {
#if defined(__SSE__)
        int data[] = {1, 2, 3, 4};
        auto v     = dsp::arrayToSIMD<int, 4>(data);
        int result[4];
        dsp::SIMDtoArray<int, 4>(result, v);
        REQUIRE(result[0] == 1);
        REQUIRE(result[1] == 2);
        REQUIRE(result[2] == 3);
        REQUIRE(result[3] == 4);

        // Aligned load and store
        int *alignedData = (int *)_mm_malloc(4 * sizeof(int), 16);
        alignedData[0]   = 1;
        alignedData[1]   = 2;
        alignedData[2]   = 3;
        alignedData[3]   = 4;
        v                = dsp::arrayToSIMD<int, 4, true>(alignedData);
        dsp::SIMDtoArray<int, 4, true>(alignedData, v);
        REQUIRE(alignedData[0] == 1);
        REQUIRE(alignedData[1] == 2);
        REQUIRE(alignedData[2] == 3);
        REQUIRE(alignedData[3] == 4);
        _mm_free(alignedData);
#endif
    }

    SECTION("SIMD<float, 4> SSE")
    {
#if defined(__SSE__)
        float data[] = {1.0f, 2.0f, 3.0f, 4.0f};
        auto v       = dsp::arrayToSIMD<float, 4>(data);
        float result[4];
        dsp::SIMDtoArray<float, 4>(result, v);
        REQUIRE(result[0] == 1.0f);
        REQUIRE(result[1] == 2.0f);
        REQUIRE(result[2] == 3.0f);
        REQUIRE(result[3] == 4.0f);

        // Aligned load and store
        auto *alignedData = (float *)_mm_malloc(4 * sizeof(float), 16);
        alignedData[0]    = 1.0f;
        alignedData[1]    = 2.0f;
        alignedData[2]    = 3.0f;
        alignedData[3]    = 4.0f;
        v                 = dsp::arrayToSIMD<float, 4, true>(alignedData);
        dsp::SIMDtoArray<float, 4, true>(alignedData, v);
        REQUIRE(alignedData[0] == 1.0f);
        REQUIRE(alignedData[1] == 2.0f);
        REQUIRE(alignedData[2] == 3.0f);
        REQUIRE(alignedData[3] == 4.0f);
        _mm_free(alignedData);
#endif
    }

    SECTION("SIMD<double, 2> SSE")
    {
#if defined(__SSE__)
        double data[] = {1.0, 2.0};
        auto v        = dsp::arrayToSIMD<double, 2>(data);
        double result[2];
        dsp::SIMDtoArray<double, 2>(result, v);
        REQUIRE(result[0] == 1.0);
        REQUIRE(result[1] == 2.0);

        // Aligned load and store
        auto *alignedData = (double *)_mm_malloc(2 * sizeof(double), 16);
        alignedData[0]    = 1.0;
        alignedData[1]    = 2.0;
        v                 = dsp::arrayToSIMD<double, 2, true>(alignedData);
        dsp::SIMDtoArray<double, 2, true>(alignedData, v);
        REQUIRE(alignedData[0] == 1.0);
        REQUIRE(alignedData[1] == 2.0);
        _mm_free(alignedData);
#endif
    }

    SECTION("SIMD<int, 8> AVX")
    {
#if defined(__AVX__)
        int data[] = {1, 2, 3, 4, 5, 6, 7, 8};
        auto v     = dsp::arrayToSIMD<int, 8>(data);
        int result[8];
        dsp::SIMDtoArray<int, 8>(result, v);
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
        aligned_data[0]   = 1;
        aligned_data[1]   = 2;
        aligned_data[2]   = 3;
        aligned_data[3]   = 4;
        aligned_data[4]   = 5;
        aligned_data[5]   = 6;
        aligned_data[6]   = 7;
        aligned_data[7]   = 8;
        v                 = dsp::arrayToSIMD<int, 8, true>(aligned_data);
        dsp::SIMDtoArray<int, 8, true>(aligned_data, v);
        REQUIRE(aligned_data[0] == 1);
        REQUIRE(aligned_data[1] == 2);
        REQUIRE(aligned_data[2] == 3);
        REQUIRE(aligned_data[3] == 4);
        REQUIRE(aligned_data[4] == 5);
        REQUIRE(aligned_data[5] == 6);
        REQUIRE(aligned_data[6] == 7);
        REQUIRE(aligned_data[7] == 8);
        _mm_free(aligned_data);
#endif
    }

    SECTION("SIMD<float, 8> AVX")
    {
#if defined(__AVX__)
        float data[] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};
        auto v       = dsp::arrayToSIMD<float, 8>(data);
        float result[8];
        dsp::SIMDtoArray<float, 8>(result, v);
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
        aligned_data[0]     = 1.0f;
        aligned_data[1]     = 2.0f;
        aligned_data[2]     = 3.0f;
        aligned_data[3]     = 4.0f;
        aligned_data[4]     = 5.0f;
        aligned_data[5]     = 6.0f;
        aligned_data[6]     = 7.0f;
        aligned_data[7]     = 8.0f;
        v                   = dsp::arrayToSIMD<float, 8, true>(aligned_data);
        dsp::SIMDtoArray<float, 8, true>(aligned_data, v);
        REQUIRE(aligned_data[0] == 1.0f);
        REQUIRE(aligned_data[1] == 2.0f);
        REQUIRE(aligned_data[2] == 3.0f);
        REQUIRE(aligned_data[3] == 4.0f);
        REQUIRE(aligned_data[4] == 5.0f);
        REQUIRE(aligned_data[5] == 6.0f);
        REQUIRE(aligned_data[6] == 7.0f);
        REQUIRE(aligned_data[7] == 8.0f);
        _mm_free(aligned_data);
#endif
    }

    SECTION("SIMD<double, 4> AVX")
    {
#if defined(__AVX__)
        double data[] = {1.0, 2.0, 3.0, 4.0};
        auto v        = dsp::arrayToSIMD<double, 4>(data);
        double result[4];
        dsp::SIMDtoArray<double, 4>(result, v);
        REQUIRE(result[0] == 1.0);
        REQUIRE(result[1] == 2.0);
        REQUIRE(result[2] == 3.0);
        REQUIRE(result[3] == 4.0);

        // Aligned load and store
        double *aligned_data = (double *)_mm_malloc(4 * sizeof(double), 32);
        aligned_data[0]      = 1.0;
        aligned_data[1]      = 2.0;
        aligned_data[2]      = 3.0;
        aligned_data[3]      = 4.0;
        v                    = dsp::arrayToSIMD<double, 4, true>(aligned_data);
        dsp::SIMDtoArray<double, 4, true>(aligned_data, v);
        REQUIRE(aligned_data[0] == 1.0);
        REQUIRE(aligned_data[1] == 2.0);
        REQUIRE(aligned_data[2] == 3.0);
        REQUIRE(aligned_data[3] == 4.0);
        _mm_free(aligned_data);
#endif
    }
}

TEST_CASE("SIMD Setting", "[simd]")
{
    SECTION("SIMD<int, 1>")
    {
        auto v = dsp::SIMD<int, 1>::set(1);
        REQUIRE(v == 1);
    }

    SECTION("SIMD<int, 4> SSE")
    {
#if defined(__SSE__)
        auto v = dsp::SIMD<int, 4>::set(1);
        int result[4];
        _mm_storeu_si128((__m128i *)result, v);
        REQUIRE(result[0] == 1);
        REQUIRE(result[1] == 1);
        REQUIRE(result[2] == 1);
        REQUIRE(result[3] == 1);
#endif
    }

    SECTION("SIMD<float, 4> SSE")
    {
#if defined(__SSE__)
        auto v = dsp::SIMD<float, 4>::set(1.0f);
        float result[4];
        _mm_storeu_ps(result, v);
        REQUIRE(result[0] == 1.0f);
        REQUIRE(result[1] == 1.0f);
        REQUIRE(result[2] == 1.0f);
        REQUIRE(result[3] == 1.0f);
#endif
    }

    SECTION("SIMD<double, 2> SSE")
    {
#if defined(__SSE__)
        auto v = dsp::SIMD<double, 2>::set(1.0);
        double result[2];
        _mm_storeu_pd(result, v);
        REQUIRE(result[0] == 1.0);
        REQUIRE(result[1] == 1.0);
#endif
    }

    SECTION("SIMD<int, 8> AVX")
    {
#if defined(__AVX__)
        auto v = dsp::SIMD<int, 8>::set(1);
        int result[8];
        _mm256_storeu_si256((__m256i *)result, v);
        REQUIRE(result[0] == 1);
        REQUIRE(result[1] == 1);
        REQUIRE(result[2] == 1);
        REQUIRE(result[3] == 1);
        REQUIRE(result[4] == 1);
        REQUIRE(result[5] == 1);
        REQUIRE(result[6] == 1);
        REQUIRE(result[7] == 1);
#endif
    }

    SECTION("SIMD<float, 8> AVX")
    {
#if defined(__AVX__)
        auto v = dsp::SIMD<float, 8>::set(1.0f);
        float result[8];
        _mm256_storeu_ps(result, v);
        REQUIRE(result[0] == 1.0f);
        REQUIRE(result[1] == 1.0f);
        REQUIRE(result[2] == 1.0f);
        REQUIRE(result[3] == 1.0f);
        REQUIRE(result[4] == 1.0f);
        REQUIRE(result[5] == 1.0f);
        REQUIRE(result[6] == 1.0f);
        REQUIRE(result[7] == 1.0f);
#endif
    }

    SECTION("SIMD<double, 4> AVX")
    {
#if defined(__AVX__)
        auto v = dsp::SIMD<double, 4>::set(1.0);
        double result[4];
        _mm256_storeu_pd(result, v);
        REQUIRE(result[0] == 1.0);
        REQUIRE(result[1] == 1.0);
        REQUIRE(result[2] == 1.0);
        REQUIRE(result[3] == 1.0);
#endif
    }
}
