#pragma once

#include <immintrin.h>

#ifndef SIMDSIZE
#if defined(__AVX2__)
#define SIMDSIZE 32
#elif defined(__SSE__)
#define SIMDSIZE 16
#else
#define SIMDSIZE 4
#endif
#endif

namespace dsp
{

// class holding the simd type, loading & storing function given Type & Size
template <typename T, int Size> struct SIMD;
// easier definition of SIMD<Type,Size>::type
template <typename T, int Size> using SIMD_t = typename SIMD<T, Size>::type;

// get type and size

// convert and array to simd given type & size
template <typename T, int Size, bool Aligned = false>
inline SIMD_t<T, Size> arrayToSIMD(const T *x)
{
    if constexpr (Aligned) return SIMD<T, Size>::load(x);
    else
        return SIMD<T, Size>::loadu(x);
}
// store a simd variable into an array
template <typename T, int Size, bool Aligned = false>
inline void SIMDtoArray(T *x, SIMD_t<T, Size> v)
{
    if constexpr (Aligned) return SIMD<T, Size>::store(x, v);
    else
        return SIMD<T, Size>::storeu(x, v);
}

/* SIMD specialization */
template <typename T> struct SIMD<T, 1> {
    using type = T;
    static constexpr auto load(const T *x) { return *x; }
    static constexpr auto loadu(const T *x) { return load(x); }
    static constexpr auto set(T x) { return x; }
    static constexpr void store(T *x, T v) { *x = v; }
    static constexpr void storeu(T *x, T v) { store(x, v); }
};

/* SSE */
#if defined(__SSE__)

/* integer */
struct SIMDint {
};

template <> struct SIMD<int, 2> {
    using type                   = __m128i;
    static constexpr auto load   = _mm_loadu_si64;
    static constexpr auto set    = _mm_set1_epi32;
    static constexpr auto loadu  = _mm_loadu_si64;
    static constexpr auto store  = _mm_storeu_si64;
    static constexpr auto storeu = _mm_storeu_si64;
};
template <> struct SIMD<int, 4> {
    using type = __m128i;
    static auto load(const int *x)
    {
        return _mm_load_si128((const __m128i *)x);
    }
    static constexpr auto set = _mm_set1_epi32;
    static auto loadu(const int *x)
    {
        return _mm_loadu_si128((const __m128i_u *)x);
    }
    static void store(int *x, type v) { _mm_store_si128((__m128i *)x, v); }
    static void storeu(int *x, type v) { _mm_storeu_si128((__m128i_u *)x, v); }
};

template <> struct SIMD<float, 2> {
    using type = float __attribute__((__vector_size__(8), __aligned__(8)));
    static constexpr type load(const float *x) { return *(const type *)x; }
    static constexpr auto loadu = load;
    static constexpr type set(float x) { return type{x, x}; }
    static void store(float *x, type v) { *(type *)x = v; }
    static constexpr auto storeu = store;
};
template <> struct SIMD<float, 4> {
    using type                   = __m128;
    static constexpr auto load   = _mm_load_ps;
    static constexpr auto loadu  = _mm_loadu_ps;
    static constexpr auto set    = _mm_set1_ps;
    static constexpr auto store  = _mm_store_ps;
    static constexpr auto storeu = _mm_storeu_ps;
};

template <> struct SIMD<double, 2> {
    using type                   = __m128d;
    static constexpr auto load   = _mm_load_pd;
    static constexpr auto loadu  = _mm_loadu_pd;
    static constexpr auto set    = _mm_set1_pd;
    static constexpr auto store  = _mm_store_pd;
    static constexpr auto storeu = _mm_storeu_pd;
};
#endif

/* AVX */
#if defined(__AVX__)
template <> struct SIMD<int, 8> {
    using type = __m256i;
    static type load(const int *x)
    {
        return _mm256_load_si256((const type *)x);
    }
    static type loadu(const int *x)
    {
        return _mm256_loadu_si256((const __m256i_u *)x);
    }
    static constexpr auto set = _mm256_set1_epi32;
    static void store(int *x, type v) { _mm256_store_si256((type *)x, v); }
    static void storeu(int *x, type v)
    {
        _mm256_storeu_si256((__m256i_u *)x, v);
    }
};
template <> struct SIMD<float, 8> {
    using type                   = __m256;
    static constexpr auto load   = _mm256_load_ps;
    static constexpr auto loadu  = _mm256_loadu_ps;
    static constexpr auto set    = _mm256_set1_ps;
    static constexpr auto store  = _mm256_store_ps;
    static constexpr auto storeu = _mm256_storeu_ps;
};

template <> struct SIMD<double, 4> {
    using type                   = __m256d;
    static constexpr auto load   = _mm256_load_pd;
    static constexpr auto loadu  = _mm256_loadu_pd;
    static constexpr auto set    = _mm256_set1_pd;
    static constexpr auto store  = _mm256_store_pd;
    static constexpr auto storeu = _mm256_storeu_pd;
};
#endif
} // namespace dsp
