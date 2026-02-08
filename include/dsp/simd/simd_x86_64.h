#pragma once

#include "../cpu/defines.h"
#include "simd_default.h"

#if DSP_X86

#include <immintrin.h>

#if DSP_AVX
#define DSP_MAX_VEC_SIZE     32
#define DSP_PADDING_VEC_SIZE 32

#elif DSP_SSE2
#define DSP_MAX_VEC_SIZE 16

// padding size should have max possible value when dispatching is active
#if DSP_X86_DISPATCH
#define DSP_PADDING_VEC_SIZE 32
#else
#define DSP_PADDING_VEC_SIZE 32
#endif

#endif

#define DSP_SIMD_DOUBLE

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T, size_t N> struct intrin;

#if DSP_SSE2

struct floatx2_t {
    floatx2_t(__m128 v) : v_(v) {}
    always_inline operator __m128() const { return v_; }

  private:
    __m128 v_;
};
struct intx2_t {
    intx2_t(__m128i v) : v_(v) {}
    always_inline operator __m128i() const { return v_; }

  private:
    __m128i v_;
};

///////////////// INT32x4 ////////////////
template <> struct intrin<int32_t, 4> {
    using type                 = __m128i;
    using basetype             = int32_t;
    using maskbasetype         = basetype;
    using masktype             = type;
    static constexpr auto init = _mm_set1_epi32;
    static auto load(const basetype *x)
    {
        return _mm_load_si128((const type *)x);
    }
    static auto loadu(const basetype *x)
    {
        return _mm_loadu_si128((const type *)x);
    }
    static void store(basetype *x, type v) { _mm_store_si128((type *)x, v); }
    static void storeu(basetype *x, type v) { _mm_storeu_si128((type *)x, v); }
    static constexpr auto add = _mm_add_epi32;
    static constexpr auto sub = _mm_sub_epi32;
    static always_inline type vectorcall neg(type x)
    {
        return sub(init(basetype(0)), x);
    }
    static always_inline type vectorcall mul(type x1, type x2)
    {
        auto mul0_2 = _mm_mul_epu32(x1, x2);
        auto mul1_3 =
            _mm_mul_epu32(_mm_shuffle_epi32(x1, _MM_SHUFFLE(2, 3, 0, 1)),
                          _mm_shuffle_epi32(x2, _MM_SHUFFLE(2, 3, 0, 1)));
        return _mm_unpacklo_epi32(
            _mm_shuffle_epi32(mul0_2, _MM_SHUFFLE(0, 0, 2, 0)),
            _mm_shuffle_epi32(mul1_3, _MM_SHUFFLE(0, 0, 2, 0)));
    }
    static constexpr auto bitAnd    = _mm_and_si128;
    static constexpr auto bitAndNot = _mm_andnot_si128;
    static constexpr auto bitOr     = _mm_or_si128;
    static constexpr auto bitXor    = _mm_xor_si128;
    static always_inline type vectorcall bitNeg(type x)
    {
        return bitAndNot(x, _mm_set_epi32(-1, -1, -1, -1));
    }

    static constexpr auto bitShiftLeft  = _mm_slli_epi32;
    static constexpr auto bitShiftRight = _mm_srai_epi32;

#if defined(__SSE4_1__)
    static constexpr auto max = _mm_max_epi32;
    static constexpr auto min = _mm_min_epi32;
#else
    static always_inline type vectorcall max(type x1, type x2)
    {
        auto mask = cmpgt(x1, x2);
        return blend(mask, x1, x2);
    }
    static always_inline type vectorcall min(type x1, type x2)
    {
        auto mask = cmpgt(x2, x1);
        return blend(mask, x1, x2);
    }
#endif

    template <size_t Id, size_t K>
    static always_inline auto vectorcall getlane(type value)
    {
        switch (Id * K) {
        default:
            break;
        case 1:
            value = _mm_shuffle_epi32(value, _MM_SHUFFLE(0, 3, 2, 1));
            break;
        case 2:
            value = _mm_shuffle_epi32(value, _MM_SHUFFLE(1, 0, 3, 2));
            break;
        case 3:
            value = _mm_shuffle_epi32(value, _MM_SHUFFLE(2, 1, 0, 3));
            break;
        }
        if constexpr (K == 1) {
            return _mm_cvtsi128_si32(value);
        } else if constexpr (K == 2) {
            return intx2_t{value};
        } else
            return value;
    }

#if defined(__SSE3__)
    static constexpr auto abs = _mm_abs_epi32;
#else
    static always_inline type vectorcall abs(type x)
    {
        auto mask = cmpgt(x, init(0));
        return blend(mask, x, neg(x));
    }
#endif

#if DSP_AVX
    template <size_t K>
    static always_inline auto vectorcall duplicate([[maybe_unused]] type value)
    {
        if constexpr (K == 2) {
            auto vlo = _mm_shuffle_epi32(value, _MM_SHUFFLE(1, 1, 0, 0));
            auto vhi = _mm_shuffle_epi32(value, _MM_SHUFFLE(3, 3, 2, 2));
            return _mm256_set_m128i(vhi, vlo);
        } else
            static_assert(false, "wrong duplicate");
    }
#endif

    template <int Shift> static always_inline type vectorcall shift(type value)
    {
        if constexpr (Shift > 0) {
            return _mm_slli_si128(value, 4 * Shift);
        } else {
            return _mm_srli_si128(value, -4 * Shift);
        }
    }

    static always_inline type vectorcall flip1(type value)
    {
        return _mm_shuffle_epi32(value, _MM_SHUFFLE(2, 3, 0, 1));
    }
    static always_inline type vectorcall flip2(type value)
    {
        return _mm_shuffle_epi32(value, _MM_SHUFFLE(1, 0, 3, 2));
    }

    static always_inline basetype vectorcall sum(type value)
    {
        type flip = flip2(value);
        type sum  = add(value, flip);
        type swap = flip1(sum);
        return _mm_cvtsi128_si32(add(sum, swap));
    }
    static always_inline intx2_t vectorcall reduce2(type value)
    {
        type flip = flip2(value);
        return add(value, flip);
    }

    static always_inline type push(type x, basetype v)
    {
        return _mm_castps_si128(
            _mm_move_ss(_mm_castsi128_ps(shift<1>(x)),
                        _mm_castsi128_ps(_mm_cvtsi32_si128(v))));
    }
    static always_inline type push(type x, intx2_t v)
    {
        return _mm_castps_si128(_mm_shuffle_ps(
            _mm_castsi128_ps(v), _mm_castsi128_ps(x), _MM_SHUFFLE(1, 0, 1, 0)));
    }

    static constexpr auto cmpeq = _mm_cmpeq_epi32;
    static constexpr auto cmpgt = _mm_cmpgt_epi32;
    static constexpr auto cmplt = _mm_cmplt_epi32;
    static always_inline masktype vectorcall cmpge(type x2, type x1)
    {
        return bitNeg(cmpgt(x1, x2));
    }
    static always_inline masktype vectorcall cmple(type x2, type x1)
    {
        return bitNeg(cmplt(x1, x2));
    }

    // blend
    static always_inline type vectorcall blend(masktype m, type x2, type x1)
    {
#if defined(__SSE4_1__)
        return _mm_blendv_epi8(x1, x2, m);
#else
        return bitOr(bitAndNot(m, x1), bitAnd(m, x2));
#endif
    }

    static constexpr auto any = _mm_movemask_epi8;
    static always_inline bool vectorcall all(masktype m)
    {
#if defined(__SSE4_1__)
        return _mm_test_all_ones(m);
#else
        return any(m) == 0xffff;
#endif
    }

    // convert
    static always_inline type convert(float value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(double value)
    {
        return init(static_cast<basetype>(value));
    }

    static always_inline type convert(intx2_t value)
    {
        return _mm_shuffle_epi32(value, _MM_SHUFFLE(1, 0, 1, 0));
    }

    static always_inline type convert(floatx2_t value)
    {
        return convert(_mm_shuffle_ps(value, value, _MM_SHUFFLE(1, 0, 1, 0)));
    }

    static always_inline type convert(__m128 value)
    {
        return _mm_cvttps_epi32(value);
    }
    static always_inline type convert(__m128d value)
    {
        auto v = _mm_cvttpd_epi32(value);
        return _mm_shuffle_epi32(v, _MM_SHUFFLE(1, 0, 1, 0));
    }

#if DSP_AVX
    static always_inline type convert(__m256i value)
    {
        return _mm256_castsi256_si128(value);
    }
    static always_inline type convert(__m256 value)
    {
        return convert(_mm256_castps256_ps128(value));
    }
    static always_inline type convert(__m256d value)
    {
        return _mm256_cvttpd_epi32(value);
    }
#endif
};

///////////////// FLOAT x 4 ////////////////
template <> struct intrin<float, 4> {
    using basetype               = float;
    using type                   = __m128;
    using maskbasetype           = basetype;
    using masktype               = type;
    static constexpr auto init   = _mm_set1_ps;
    static constexpr auto load   = _mm_load_ps;
    static constexpr auto loadu  = _mm_loadu_ps;
    static constexpr auto store  = _mm_store_ps;
    static constexpr auto storeu = _mm_storeu_ps;
    static constexpr auto add    = _mm_add_ps;
    static constexpr auto sub    = _mm_sub_ps;
    static always_inline type vectorcall neg(type x)
    {
        return sub(init(basetype(0)), x);
    }
    static constexpr auto mul  = _mm_mul_ps;
    static constexpr auto div  = _mm_div_ps;
    static constexpr auto sqrt = _mm_sqrt_ps;

    static always_inline masktype vectorcall signbit(type x)
    {
        return (masktype)_mm_srai_epi32((__m128i)x, 31);
    }

    static constexpr auto bitAnd    = _mm_and_ps;
    static constexpr auto bitAndNot = _mm_andnot_ps;
    static constexpr auto bitOr     = _mm_or_ps;
    static constexpr auto bitXor    = _mm_xor_ps;
    static always_inline type vectorcall bitNeg(type x)
    {
        return bitAndNot(x, (__m128)_mm_set_epi32(-1, -1, -1, -1));
    }
    static constexpr auto max = _mm_max_ps;
    static constexpr auto min = _mm_min_ps;

    template <size_t Id, size_t K>
    static always_inline auto vectorcall getlane(type value)
    {
        switch (Id * K) {
        default:
            break;
        case 1:
            value = _mm_shuffle_ps(value, value, _MM_SHUFFLE(0, 3, 2, 1));
            break;
        case 2:
            value = _mm_shuffle_ps(value, value, _MM_SHUFFLE(1, 0, 3, 2));
            break;
        case 3:
            value = _mm_shuffle_ps(value, value, _MM_SHUFFLE(2, 1, 0, 3));
            break;
        }
        if constexpr (K == 1) {
            return _mm_cvtss_f32(value);
        } else if constexpr (K == 2) {
            return floatx2_t{value};
        } else
            return value;
    }

    static always_inline type vectorcall abs(type value)
    {
        return bitAnd(value, init(floatMask<basetype>::kNotSign.f));
    }

    template <int Shift> static always_inline type vectorcall shift(type value)
    {
        __m128i r;
        if constexpr (Shift > 0) {
            r = _mm_slli_si128(_mm_castps_si128(value), 4 * Shift);
        } else {
            r = _mm_srli_si128(_mm_castps_si128(value), -4 * Shift);
        }
        return _mm_castsi128_ps(r);
    }

#if DSP_AVX
    template <size_t K>
    static always_inline auto vectorcall duplicate([[maybe_unused]] type value)
    {
        if constexpr (K == 2) {
            auto vlo = _mm_shuffle_ps(value, value, _MM_SHUFFLE(1, 1, 0, 0));
            auto vhi = _mm_shuffle_ps(value, value, _MM_SHUFFLE(3, 3, 2, 2));
            return _mm256_set_m128(vhi, vlo);
        } else
            static_assert(false, "wrong duplicate");
    }
#endif

    static always_inline type vectorcall flip1(type value)
    {
        return _mm_shuffle_ps(value, value, _MM_SHUFFLE(2, 3, 0, 1));
    }
    static always_inline type vectorcall flip2(type value)
    {
        return _mm_shuffle_ps(value, value, _MM_SHUFFLE(1, 0, 3, 2));
    }
    static always_inline basetype vectorcall sum(type value)
    {
        type flip = flip2(value);
        type sum  = add(value, flip);
        type swap = flip1(sum);
        return _mm_cvtss_f32(add(sum, swap));
    }
    static always_inline floatx2_t vectorcall reduce2(type value)
    {
        type flip = flip2(value);
        return add(value, flip);
    }

    static always_inline type push(type x, basetype v)
    {
        return _mm_move_ss(shift<1>(x), _mm_set_ss(v));
    }
    static always_inline type push(type x, floatx2_t v)
    {
        return _mm_shuffle_ps(v, x, _MM_SHUFFLE(1, 0, 1, 0));
    }

    static constexpr auto cmpeq = _mm_cmpeq_ps;
    static constexpr auto cmpgt = _mm_cmpgt_ps;
    static constexpr auto cmplt = _mm_cmplt_ps;
    static constexpr auto cmpge = _mm_cmpge_ps;
    static constexpr auto cmple = _mm_cmple_ps;

    static always_inline type vectorcall blend(masktype m, type x2, type x1)
    {
#if defined(__SSE4_1__)
        return _mm_blendv_ps(x1, x2, m);
#else
        return bitOr(bitAndNot(m, x1), bitAnd(m, x2));
#endif
    }

    static constexpr auto any = _mm_movemask_ps;
    static always_inline bool vectorcall all(masktype m)
    {
#if defined(__SSE4_1__)
        return _mm_test_all_ones((__m128i)m);
#else
        return any(m) == 0x000f;
#endif
    }

    // convert
    static always_inline type convert(float value) { return init(value); }
    static always_inline type convert(floatx2_t value)
    {
        return _mm_shuffle_ps(value, value, _MM_SHUFFLE(1, 0, 1, 0));
    }
    static always_inline type convert(intx2_t value)
    {
        return convert(_mm_shuffle_epi32(value, _MM_SHUFFLE(1, 0, 1, 0)));
    }
    static always_inline type convert(__m128i value)
    {
        return _mm_cvtepi32_ps(value);
    }
    static always_inline type convert(__m128d value)
    {
        return convert(static_cast<floatx2_t>(_mm_cvtpd_ps(value)));
    }

#if DSP_AVX
    static always_inline type convert(__m256 value)
    {
        return _mm256_castps256_ps128(value);
    }
    static always_inline type convert(__m256i value)
    {
        return convert(intrin<int, 4>::convert(value));
    }
    static always_inline type convert(__m256d value)
    {
        return _mm256_cvtpd_ps(value);
    }
#endif
};

///////////////// FLOAT x 2 ////////////////
template <> struct intrin<float, 2> {
    using base         = intrin<float, 4>;
    using basetype     = float;
    using type         = floatx2_t;
    using maskbasetype = basetype;
    using masktype     = type;

    // static constexpr auto init = base::init;
    static always_inline type vectorcall init(basetype v)
    {
        return _mm_set_ps(0, 0, v, v);
    }
    static always_inline type vectorcall loadu(const basetype *src)
    {
        return ((__m128)_mm_loadl_epi64((const __m128i *)src));
    }
    static constexpr auto load = loadu;
    static always_inline void vectorcall storeu(basetype *dest, type value)
    {
        _mm_storel_epi64((__m128i *)dest, (__m128i) static_cast<__m128>(value));
    }
    static constexpr auto store = storeu;
    static always_inline type vectorcall add(type x1, type x2)
    {
        return base::add(x1, x2);
    }
    static always_inline type vectorcall sub(type x1, type x2)
    {
        return base::sub(x1, x2);
    }
    static always_inline type vectorcall neg(type x) { return base::neg(x); }
    static always_inline type vectorcall mul(type x1, type x2)
    {
        return base::mul(x1, x2);
    }
    static always_inline type vectorcall div(type x1, type x2)
    {
        return base::div(x1, x2);
    }
    static always_inline type vectorcall sqrt(type x) { return base::sqrt(x); }
    static always_inline type vectorcall signbit(type x)
    {
        return base::signbit(x);
    }
    static always_inline type vectorcall bitAnd(type x1, type x2)
    {
        return base::bitAnd(x1, x2);
    }
    static always_inline type vectorcall bitOr(type x1, type x2)
    {
        return base::bitOr(x1, x2);
    }
    static always_inline type vectorcall bitXor(type x1, type x2)
    {
        return base::bitXor(x1, x2);
    }
    static always_inline type vectorcall bitNeg(type x)
    {
        return base::bitNeg(x);
    }
    static always_inline type vectorcall max(type x1, type x2)
    {
        return base::max(x1, x2);
    }
    static always_inline type vectorcall min(type x1, type x2)
    {
        return base::min(x1, x2);
    }

    template <size_t Id, size_t K>
    static always_inline auto vectorcall getlane(type value)
    {
        if (Id == 1) {
            value = _mm_shuffle_ps(value, value, _MM_SHUFFLE(0, 3, 2, 1));
        }
        if constexpr (K == 1) {
            return _mm_cvtss_f32(value);
        } else
            return value;
    }

    static always_inline type vectorcall abs(type x) { return base::abs(x); }

    template <int Shift> static always_inline type vectorcall shift(type value)
    {
        if constexpr (Shift > 0) {
            return _mm_shuffle_ps(value, value, _MM_SHUFFLE(3, 2, 0, 3));
        } else {
            return _mm_shuffle_ps(value, value, _MM_SHUFFLE(3, 2, 3, 1));
        }
    }

    template <size_t K>
    static always_inline auto vectorcall duplicate([[maybe_unused]] type value)
    {
        if constexpr (K == 2) {
            return _mm_shuffle_ps(value, value, _MM_SHUFFLE(1, 1, 0, 0));
        }
#if DSP_AVX
        else if constexpr (K == 4) {
            auto vlo = _mm_shuffle_ps(value, value, _MM_SHUFFLE(0, 0, 0, 0));
            auto vhi = _mm_shuffle_ps(value, value, _MM_SHUFFLE(1, 1, 1, 1));
            return _mm256_set_m128(vhi, vlo);
        }
#endif
        else
            static_assert(false, "wrong duplicate");
    }

    static always_inline type vectorcall flip1(type value)
    {
        return base::flip1(value);
    }
    static always_inline basetype vectorcall sum(type value)
    {
        type swap = flip1(value);
        type sum  = add(value, swap);
        return _mm_cvtss_f32(sum);
    }

    static always_inline type push(type x, basetype v)
    {
        return _mm_move_ss(shift<1>(x), _mm_set_ss(v));
    }

    static always_inline type vectorcall cmpeq(type x1, type x2)
    {
        return base::cmpeq(x1, x2);
    }
    static always_inline type vectorcall cmpgt(type x1, type x2)
    {
        return base::cmpgt(x1, x2);
    }
    static always_inline type vectorcall cmplt(type x1, type x2)
    {
        return base::cmplt(x1, x2);
    }
    static always_inline type vectorcall cmpge(type x1, type x2)
    {
        return base::cmpge(x1, x2);
    }
    static always_inline type vectorcall cmple(type x1, type x2)
    {
        return base::cmple(x1, x2);
    }
    static always_inline type vectorcall blend(masktype m, type x1, type x2)
    {
        return base::blend(m, x1, x2);
    }
    static always_inline auto vectorcall any(masktype m)
    {
        return base::any(m);
    }
    static always_inline auto vectorcall all(masktype m)
    {
        return (any(m) & 0x0003) == 0x0003;
    }

    // convert
    static always_inline type vectorcall convert(float x) { return init(x); }
    static always_inline type convert(intx2_t value)
    {
        return _mm_cvtepi32_ps(value);
    }
    static always_inline type convert(__m128d value)
    {
        return _mm_cvtpd_ps(value);
    }
};

///////////////// INT x 2 ////////////////
template <> struct intrin<int, 2> {
    using base         = intrin<int, 4>;
    using basetype     = int;
    using type         = intx2_t;
    using maskbasetype = basetype;
    using masktype     = type;

    static constexpr auto init = base::init;
    static always_inline type vectorcall loadu(const basetype *src)
    {
        return ((__m128i)_mm_load_sd((const double *)src));
    }
    static constexpr auto load = loadu;
    static always_inline void vectorcall storeu(basetype *dest, type value)
    {
        _mm_store_sd((double *)dest, (__m128d) static_cast<__m128i>(value));
    }
    static constexpr auto store = storeu;
    static always_inline type vectorcall add(type x1, type x2)
    {
        return base::add(x1, x2);
    }
    static always_inline type vectorcall sub(type x1, type x2)
    {
        return base::sub(x1, x2);
    }
    static always_inline type vectorcall neg(type x) { return base::neg(x); }
    static always_inline type vectorcall mul(type x1, type x2)
    {
        return base::mul(x1, x2);
    }
    static always_inline type vectorcall bitAnd(type x1, type x2)
    {
        return base::bitAnd(x1, x2);
    }
    static always_inline type vectorcall bitOr(type x1, type x2)
    {
        return base::bitOr(x1, x2);
    }
    static always_inline type vectorcall bitXor(type x1, type x2)
    {
        return base::bitXor(x1, x2);
    }
    static always_inline type vectorcall bitNeg(type x)
    {
        return base::bitNeg(x);
    }
    static always_inline type vectorcall max(type x1, type x2)
    {
        return base::max(x1, x2);
    }

    static always_inline type vectorcall bitShiftLeft(type x1, int shift)
    {
        return base::bitShiftLeft(x1, shift);
    }
    static always_inline type vectorcall bitShiftRight(type x1, int shift)
    {
        return base::bitShiftRight(x1, shift);
    }

    static always_inline type vectorcall min(type x1, type x2)
    {
        return base::min(x1, x2);
    }

    template <size_t Id, size_t K>
    static always_inline auto vectorcall getlane(type value)
    {
        if (Id == 1) {
            value = _mm_shuffle_epi32(value, _MM_SHUFFLE(0, 3, 2, 1));
        }
        if constexpr (K == 1) {
            return _mm_cvtsi128_si32(value);
        } else
            return value;
    }

    static always_inline type vectorcall abs(type x) { return base::abs(x); }

    template <int Shift> static always_inline type vectorcall shift(type value)
    {
        if constexpr (Shift > 0) {
            return _mm_shuffle_epi32(value, _MM_SHUFFLE(3, 2, 0, 3));
        } else {
            return _mm_shuffle_epi32(value, _MM_SHUFFLE(3, 2, 3, 1));
        }
    }

    template <size_t K>
    static always_inline auto vectorcall duplicate([[maybe_unused]] type value)
    {
        if constexpr (K == 2) {
            return _mm_shuffle_epi32(value, _MM_SHUFFLE(1, 1, 0, 0));
        }
#if DSP_AVX
        else if constexpr (K == 4) {
            auto vlo = _mm_shuffle_epi32(value, _MM_SHUFFLE(0, 0, 0, 0));
            auto vhi = _mm_shuffle_epi32(value, _MM_SHUFFLE(1, 1, 1, 1));
            return _mm256_set_m128i(vhi, vlo);
        }
#endif
        else
            static_assert(false, "wrong duplicate");
    }

    static always_inline type vectorcall flip1(type value)
    {
        return base::flip1(value);
    }
    static always_inline basetype vectorcall sum(type value)
    {
        type swap = flip1(value);
        type sum  = add(value, swap);
        return _mm_cvtsi128_si32(sum);
    }

    static always_inline type push(type x, basetype v)
    {
        return _mm_castps_si128(
            _mm_move_ss(_mm_castsi128_ps(shift<1>(x)),
                        _mm_castsi128_ps(_mm_cvtsi32_si128(v))));
    }

    static always_inline type vectorcall cmpeq(type x1, type x2)
    {
        return base::cmpeq(x1, x2);
    }
    static always_inline type vectorcall cmpgt(type x1, type x2)
    {
        return base::cmpgt(x1, x2);
    }
    static always_inline type vectorcall cmplt(type x1, type x2)
    {
        return base::cmplt(x1, x2);
    }
    static always_inline type vectorcall cmpge(type x1, type x2)
    {
        return base::cmpge(x1, x2);
    }
    static always_inline type vectorcall cmple(type x1, type x2)
    {
        return base::cmple(x1, x2);
    }
    static always_inline type vectorcall blend(masktype m, type x1, type x2)
    {
        return base::blend(m, x1, x2);
    }
    static always_inline auto vectorcall any(masktype m)
    {
        return base::any(m) & 0x00ff;
    }
    static always_inline auto vectorcall all(masktype m)
    {
        return any(m) == 0x00ff;
    }

    // convert
    static always_inline type vectorcall convert(float x) { return init(x); }
    static always_inline type vectorcall convert(floatx2_t x)
    {
        return _mm_cvttps_epi32(x);
    }
    static always_inline type convert(__m128d value)
    {
        return _mm_cvttpd_epi32(value);
    }
};

///////////////// DOUBLE x 2 ////////////////
template <> struct intrin<double, 2> {
    using basetype               = double;
    using type                   = __m128d;
    using maskbasetype           = basetype;
    using masktype               = type;
    static constexpr auto init   = _mm_set1_pd;
    static constexpr auto load   = _mm_load_pd;
    static constexpr auto loadu  = _mm_loadu_pd;
    static constexpr auto store  = _mm_store_pd;
    static constexpr auto storeu = _mm_storeu_pd;
    static constexpr auto add    = _mm_add_pd;
    static constexpr auto sub    = _mm_sub_pd;
    static always_inline type vectorcall neg(type x)
    {
        return sub(init(basetype(0)), x);
    }
    static constexpr auto mul  = _mm_mul_pd;
    static constexpr auto div  = _mm_div_pd;
    static constexpr auto sqrt = _mm_sqrt_pd;

    static always_inline masktype vectorcall signbit(type x)
    {
        auto shift = _mm_srai_epi32((__m128i)x, 31);
        return (masktype)_mm_shuffle_epi32(shift, _MM_SHUFFLE(3, 3, 1, 1));
    }

    static constexpr auto bitAnd    = _mm_and_pd;
    static constexpr auto bitAndNot = _mm_andnot_pd;
    static constexpr auto bitOr     = _mm_or_pd;
    static constexpr auto bitXor    = _mm_xor_pd;
    static always_inline type vectorcall bitNeg(type x)
    {
        return bitAndNot(x, (__m128d)_mm_set_epi32(-1, -1, -1, -1));
    }
    static constexpr auto max = _mm_max_pd;
    static constexpr auto min = _mm_min_pd;

    template <size_t Id, size_t K>
    static always_inline auto vectorcall getlane(type value)
    {
        if (Id == 1) {
            value = _mm_shuffle_pd(value, value, _MM_SHUFFLE2(0, 1));
        }
        if constexpr (K == 1) {
            return _mm_cvtsd_f64(value);
        } else
            return value;
    }

    static always_inline type vectorcall abs(type value)
    {
        return bitAnd(value, init(*reinterpret_cast<const basetype *>(
                                 &floatMask<basetype>::kNotSign)));
    }

#if DSP_AVX
    template <size_t K>
    static always_inline auto vectorcall duplicate([[maybe_unused]] type value)
    {
        if constexpr (K == 2) {
            auto vlo = _mm_shuffle_pd(value, value, _MM_SHUFFLE2(0, 0));
            auto vhi = _mm_shuffle_pd(value, value, _MM_SHUFFLE2(1, 1));
            return _mm256_set_m128d(vhi, vlo);
        } else
            static_assert(false, "wrong duplicate");
    }
#endif

    template <int Shift> static always_inline type vectorcall shift(type value)
    {
        __m128i r;
        if constexpr (Shift > 0) {
            r = _mm_slli_si128(_mm_castpd_si128(value), 8 * Shift);
        } else {
            r = _mm_srli_si128(_mm_castpd_si128(value), -8 * Shift);
        }
        return _mm_castsi128_pd(r);
    }

    static always_inline type vectorcall flip1(type value)
    {
        return _mm_shuffle_pd(value, value, _MM_SHUFFLE2(0, 1));
    }
    static always_inline basetype vectorcall sum(type value)
    {
        type flip = flip1(value);
        return _mm_cvtsd_f64(add(value, flip));
    }

    static always_inline type push(type x, basetype v)
    {
        return _mm_move_sd(shift<1>(x), _mm_set_sd(v));
    }

    static constexpr auto cmpeq = _mm_cmpeq_pd;
    static constexpr auto cmpgt = _mm_cmpgt_pd;
    static constexpr auto cmplt = _mm_cmplt_pd;
    static constexpr auto cmpge = _mm_cmpge_pd;
    static constexpr auto cmple = _mm_cmple_pd;

    static always_inline type vectorcall blend(masktype m, type x2, type x1)
    {
#if defined(__SSE4_1__)
        return _mm_blendv_pd(x1, x2, m);
#else
        return bitOr(bitAndNot(m, x1), bitAnd(m, x2));
#endif
    }

    static constexpr auto any = _mm_movemask_pd;
    static always_inline bool vectorcall all(masktype m)
    {
#if defined(__SSE4_1__)
        return _mm_test_all_ones((__m128i)m);
#else
        return any(m) == 0x0003;
#endif
    }

    // convert
    static always_inline type convert(float value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(int value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(double value) { return init(value); }

    static always_inline type convert(floatx2_t value)
    {
        return convert(static_cast<__m128>(value));
    }

    static always_inline type convert(intx2_t value)
    {
        return convert(static_cast<__m128i>(value));
    }

    static always_inline type convert(__m128i value)
    {
        return _mm_cvtepi32_pd(value);
    }
    static always_inline type convert(__m128 value)
    {
        return _mm_cvtps_pd(value);
    }

#if DSP_AVX
    static always_inline type convert(__m256d value)
    {
        return _mm256_castpd256_pd128(value);
    }
    static always_inline type convert(__m256i value)
    {
        return convert(intrin<int, 4>::convert(value));
    }
    static always_inline type convert(__m256 value)
    {
        return convert(intrin<float, 4>::convert(value));
    }
#endif
};
#endif

#if DSP_AVX
///////////////// INT32x8 ////////////////
template <> struct intrin<int32_t, 8> {
    using type                 = __m256i;
    using basetype             = int32_t;
    using maskbasetype         = basetype;
    using masktype             = type;
    static constexpr auto init = _mm256_set1_epi32;
    static auto load(const basetype *x)
    {
        return _mm256_load_si256((const type *)x);
    }
    static auto loadu(const basetype *x)
    {
        return _mm256_loadu_si256((const type *)x);
    }
    static void store(basetype *x, type v) { _mm256_store_si256((type *)x, v); }
    static void storeu(basetype *x, type v)
    {
        _mm256_storeu_si256((type *)x, v);
    }
    static constexpr auto add = _mm256_add_epi32;
    static constexpr auto sub = _mm256_sub_epi32;
    static always_inline type vectorcall neg(type x)
    {
        return sub(init(basetype(0)), x);
    }
    static always_inline type vectorcall mul(type x1, type x2)
    {
        auto mul0_2 = _mm256_mul_epi32(x1, x2);
        auto mul1_3 =
            _mm256_mul_epi32(_mm256_shuffle_epi32(x1, _MM_SHUFFLE(2, 3, 0, 1)),
                             _mm256_shuffle_epi32(x2, _MM_SHUFFLE(2, 3, 0, 1)));
        return _mm256_unpacklo_epi32(
            _mm256_shuffle_epi32(mul0_2, _MM_SHUFFLE(0, 0, 2, 0)),
            _mm256_shuffle_epi32(mul1_3, _MM_SHUFFLE(0, 0, 2, 0)));
    }
    static constexpr auto bitAnd = _mm256_and_si256;
    static constexpr auto bitOr  = _mm256_or_si256;
    static constexpr auto bitXor = _mm256_xor_si256;
    static always_inline type vectorcall bitNeg(type x)
    {
        return bitXor(init(basetype(-1)), x);
    }
    static constexpr auto bitShiftLeft  = _mm256_slli_epi32;
    static constexpr auto bitShiftRight = _mm256_srai_epi32;

    static constexpr auto max = _mm256_max_epi32;
    static constexpr auto min = _mm256_min_epi32;
    static constexpr auto abs = _mm256_abs_epi32;

    template <size_t Id, size_t K>
    static always_inline auto vectorcall getlane(type value)
    {
        switch ((Id * K) % 4) {
        default:
            break;
        case 1:
            value = _mm256_shuffle_epi32(value, _MM_SHUFFLE(0, 3, 2, 1));
            break;
        case 2:
            value = _mm256_shuffle_epi32(value, _MM_SHUFFLE(1, 0, 3, 2));
            break;
        case 3:
            value = _mm256_shuffle_epi32(value, _MM_SHUFFLE(2, 1, 0, 3));
            break;
        }

        auto r = _mm256_extractf128_si256(value, Id * K >= 4 ? 1 : 0);
        if constexpr (K == 1) {
            return _mm_cvtsi128_si32(r);
        } else if constexpr (K == 2) {
            return intx2_t{r};
        } else if constexpr (K == 4) {
            return r;
        } else
            return value;
    }

    template <int Shift> static always_inline type vectorcall shift(type value)
    {
        auto vp = _mm256_permute2f128_si256(
            value, value, Shift > 0 ? 1 << 3 : (1) + (1 << 7));

        constexpr auto kShift = Shift > 0 ? Shift : -Shift;
        if constexpr (kShift == 4) return vp;

        constexpr __m256i kZero{};
        constexpr auto kShiftAlign = kShift < 4 ? kShift : kShift - 4;

        auto v1 = Shift > 0 ? value : vp;
        auto v2 = Shift > 0 ? vp : value;

        auto va1 = kShift < 4 ? v1 : v2;
        auto va2 = kShift < 4 ? v2 : kZero;

        constexpr auto kAlign = Shift > 0 ? (kShiftAlign == 1   ? 12
                                             : kShiftAlign == 2 ? 8
                                             : kShiftAlign == 3 ? 4
                                                                : 0)
                                          : kShift * 4;
        return Shift < -4 ? _mm256_alignr_epi8(vp, kZero, kAlign)
                          : _mm256_alignr_epi8(va1, va2, kAlign);
    }

    static always_inline type vectorcall flip1(type value)
    {
        return _mm256_shuffle_epi32(value, _MM_SHUFFLE(2, 3, 0, 1));
    }
    static always_inline type vectorcall flip2(type value)
    {
        return _mm256_shuffle_epi32(value, _MM_SHUFFLE(1, 0, 3, 2));
    }
    static always_inline type vectorcall flip4(type value)
    {
        return _mm256_permute2f128_si256(value, value, _MM_SHUFFLE2(0, 3));
    }
    static always_inline basetype vectorcall sum(type value)
    {
        type flip = flip2(value);
        type sum  = add(value, flip);
        flip      = flip1(sum);
        sum       = add(sum, flip);
        flip      = flip4(sum);
        sum       = add(sum, flip);
        return _mm256_extract_epi32(sum, 0);
    }
    static always_inline intx2_t vectorcall reduce2(type value)
    {
        type flip = flip2(value);
        type sum  = add(value, flip);
        flip      = flip4(sum);
        sum       = add(sum, flip);
        return _mm256_extracti128_si256(sum, 0);
    }
    static always_inline __m128i vectorcall reduce4(type value)
    {
        auto flip = flip4(value);
        auto sum  = add(value, flip);
        return _mm256_extracti128_si256(sum, 0);
    }

    template <size_t K> static always_inline type vectorcall prefix(type value)
    {
        static_assert(K == 1 || K == 2 || K == 4);
        if (K == 1) {
            auto shift1 = _mm256_slli_si256(value, 4);
            value       = _mm256_add_epi32(value, shift1);
        }
        if (K <= 2) {
            auto shift2 = _mm256_slli_si256(value, 8);
            value       = _mm256_add_epi32(value, shift2);
        }

        __m128i sum1block;
        switch (K) {
        default:
        case 1:
            sum1block = _mm_set1_epi32(_mm256_extract_epi32(value, 3));
            break;
        case 2:
            sum1block = _mm_shuffle_epi32(_mm256_castsi256_si128(value),
                                          _MM_SHUFFLE(3, 2, 3, 2));
            break;
        case 4:
            sum1block = _mm256_castsi256_si128(value);
            break;
        }

        return _mm256_add_epi32(value, _mm256_setr_m128i(__m128i{}, sum1block));
    }

    static always_inline type push(type x, basetype v)
    {
        auto sh = _mm256_castsi256_ps(shift<1>(x));
        auto lo = _mm_move_ss(_mm256_castps256_ps128(sh),
                              _mm_castsi128_ps(_mm_cvtsi32_si128(v)));
        return _mm256_castps_si256(_mm256_insertf128_ps(sh, lo, 0));
    }
    static always_inline type push(type x, intx2_t v)
    {
        auto sh = _mm256_castsi256_ps(shift<2>(x));
        auto lo =
            _mm_shuffle_ps(_mm_castsi128_ps(v), _mm256_castps256_ps128(sh),
                           _MM_SHUFFLE(3, 2, 1, 0));
        return _mm256_castps_si256(_mm256_insertf128_ps(sh, lo, 0));
    }
    static always_inline type push(type x, __m128i v)
    {
        auto sh = shift<4>(x);
        return _mm256_inserti128_si256(sh, v, 0);
    }

    static constexpr auto cmpeq = _mm256_cmpeq_epi32;
    static constexpr auto cmpgt = _mm256_cmpgt_epi32;
    static always_inline masktype vectorcall cmpge(type x2, type x1)
    {
        // x2 >= x1 <==> ~(x1 > x2)
        return bitNeg(cmpgt(x1, x2));
    }
    static always_inline masktype vectorcall cmplt(type x2, type x1)
    {
        return cmpgt(x1, x2);
    }
    static always_inline masktype vectorcall cmple(type x2, type x1)
    {
        return cmpge(x1, x2);
    }

    static always_inline type vectorcall blend(masktype m, type x2, type x1)
    {
        return _mm256_blendv_epi8(x1, x2, m);
    }

    static constexpr auto any = _mm256_movemask_epi8;
    static always_inline bool vectorcall all(masktype m)
    {
        return any(m) == -1;
    }

    // convert
    static always_inline type convert(basetype value) { return init(value); }
    static always_inline type convert(float value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(double value)
    {
        return init(static_cast<basetype>(value));
    }

    static always_inline type convert(floatx2_t value)
    {
        return convert(intrin<float, 4>::convert(value));
    }

    static always_inline type convert(__m128i value)
    {
        return _mm256_set_m128i(value, value);
    }
    static always_inline type convert(__m128 value)
    {
        return convert(intrin<int, 4>::convert(value));
    }
    static always_inline type convert(__m128d value)
    {
        return convert(intrin<int, 4>::convert(value));
    }

    static always_inline type convert(__m256 value)
    {
        return _mm256_cvttps_epi32(value);
    }
    static always_inline type convert(__m256d value)
    {
        return convert(_mm256_cvttpd_epi32(value));
    }
};

///////////////// FLOAT x 8 ////////////////
template <> struct intrin<float, 8> {
    using basetype               = float;
    using type                   = __m256;
    using maskbasetype           = basetype;
    using masktype               = type;
    static constexpr auto init   = _mm256_set1_ps;
    static constexpr auto load   = _mm256_load_ps;
    static constexpr auto loadu  = _mm256_loadu_ps;
    static constexpr auto store  = _mm256_store_ps;
    static constexpr auto storeu = _mm256_storeu_ps;
    static constexpr auto add    = _mm256_add_ps;
    static constexpr auto sub    = _mm256_sub_ps;
    static always_inline type vectorcall neg(type x)
    {
        return sub(init(basetype(0)), x);
    }
    static constexpr auto mul  = _mm256_mul_ps;
    static constexpr auto div  = _mm256_div_ps;
    static constexpr auto sqrt = _mm256_sqrt_ps;

    static always_inline masktype vectorcall signbit(type x)
    {
        return (masktype)_mm256_srai_epi32((__m256i)x, 31);
    }

    static constexpr auto bitAnd = _mm256_and_ps;
    static constexpr auto bitOr  = _mm256_or_ps;
    static constexpr auto bitXor = _mm256_xor_ps;
    static always_inline type vectorcall bitNeg(type x)
    {
        return bitXor(init(basetype(0)), x);
    }
    static constexpr auto max = _mm256_max_ps;
    static constexpr auto min = _mm256_min_ps;

    static always_inline type vectorcall abs(type value)
    {
        return bitAnd(value, init(*reinterpret_cast<const basetype *>(
                                 &floatMask<basetype>::kNotSign)));
    }

    template <size_t Id, size_t K>
    static always_inline auto vectorcall getlane(type value)
    {
        switch ((Id * K) % 4) {
        default:
            break;
        case 1:
            value = _mm256_shuffle_ps(value, value, _MM_SHUFFLE(0, 3, 2, 1));
            break;
        case 2:
            value = _mm256_shuffle_ps(value, value, _MM_SHUFFLE(1, 0, 3, 2));
            break;
        case 3:
            value = _mm256_shuffle_ps(value, value, _MM_SHUFFLE(2, 1, 0, 3));
            break;
        }

        auto r = _mm256_extractf128_ps(value, Id * K >= 4 ? 1 : 0);
        if constexpr (K == 1) {
            return _mm_cvtss_f32(r);
        } else if constexpr (K == 2) {
            return floatx2_t{r};
        } else if constexpr (K == 4) {
            return r;
        } else
            return value;
    }

    template <int Shift> static always_inline type vectorcall shift(type value)
    {
        return _mm256_castsi256_ps(
            intrin<int, 8>::template shift<Shift>(_mm256_castps_si256(value)));
    }

    static always_inline type vectorcall flip1(type value)
    {
        return _mm256_shuffle_ps(value, value, _MM_SHUFFLE(2, 3, 0, 1));
    }
    static always_inline type vectorcall flip2(type value)
    {
        return _mm256_shuffle_ps(value, value, _MM_SHUFFLE(1, 0, 3, 2));
    }
    static always_inline type vectorcall flip4(type value)
    {
        return _mm256_permute2f128_ps(value, value, _MM_SHUFFLE2(0, 3));
    }
    static always_inline basetype vectorcall sum(type value)
    {
        type flip = flip2(value);
        type sum  = add(value, flip);
        flip      = flip1(sum);
        sum       = add(sum, flip);
        flip      = flip4(sum);
        sum       = add(sum, flip);
        return _mm256_cvtss_f32(sum);
    }

    static always_inline floatx2_t vectorcall reduce2(type value)
    {
        type flip = flip2(value);
        type sum  = add(value, flip);
        flip      = flip4(sum);
        sum       = add(sum, flip);
        return _mm256_extractf128_ps(sum, 0);
    }

    static always_inline __m128 vectorcall reduce4(type value)
    {
        auto flip = flip4(value);
        auto sum  = add(value, flip);
        return _mm256_extractf128_ps(sum, 0);
    }

    template <size_t K> static always_inline type vectorcall prefix(type value)
    {
        static_assert(K == 1 || K == 2 || K == 4);
        if (K == 1) {
            auto shift1 = _mm256_castsi256_ps(
                _mm256_slli_si256(_mm256_castps_si256(value), 4));
            value = _mm256_add_ps(value, shift1);
        }
        if (K <= 2) {
            auto shift2 = _mm256_castsi256_ps(
                _mm256_slli_si256(_mm256_castps_si256(value), 8));
            value = _mm256_add_ps(value, shift2);
        }

        __m128 sum1block;
        switch (K) {
        default:
        case 1:
            sum1block = _mm_set1_ps(value[3]);
            break;
        case 2:
            sum1block = _mm_shuffle_ps(_mm256_castps256_ps128(value),
                                       _mm256_castps256_ps128(value),
                                       _MM_SHUFFLE(3, 2, 3, 2));
            break;
        case 4:
            sum1block = _mm256_castps256_ps128(value);
            break;
        }

        return _mm256_add_ps(value, _mm256_setr_m128(__m128{}, sum1block));
    }

    static always_inline type push(type x, basetype v)
    {
        auto sh = shift<1>(x);
        auto lo = _mm_move_ss(_mm256_castps256_ps128(sh), _mm_set_ss(v));
        return _mm256_insertf128_ps(sh, lo, 0);
    }
    static always_inline type push(type x, floatx2_t v)
    {
        auto sh = shift<2>(x);
        auto lo = _mm_shuffle_ps(v, _mm256_castps256_ps128(sh),
                                 _MM_SHUFFLE(3, 2, 1, 0));
        return _mm256_insertf128_ps(sh, lo, 0);
    }
    static always_inline type push(type x, __m128 v)
    {
        auto sh = shift<4>(x);
        return _mm256_insertf128_ps(sh, v, 0);
    }

    static always_inline type vectorcall cmpeq(type x1, type x2)
    {
        return _mm256_cmp_ps(x1, x2, 0x00);
    }
    static always_inline type vectorcall cmpgt(type x1, type x2)
    {
        return _mm256_cmp_ps(x1, x2, 0x0E);
    }
    static always_inline type vectorcall cmpge(type x1, type x2)
    {
        return _mm256_cmp_ps(x1, x2, 0x0D);
    }
    static always_inline type vectorcall cmplt(type x1, type x2)
    {
        return _mm256_cmp_ps(x1, x2, 0x11);
    }
    static always_inline type vectorcall cmple(type x1, type x2)
    {
        return _mm256_cmp_ps(x1, x2, 0x12);
    }

    static always_inline type vectorcall blend(masktype m, type x2, type x1)
    {
        return _mm256_blendv_ps(x1, x2, m);
    }

    static constexpr auto any = _mm256_movemask_ps;
    static always_inline bool vectorcall all(masktype m)
    {
        return any(m) == 0x00ff;
    }

    // convert
    static always_inline type convert(float value) { return init(value); }
    static always_inline type convert(double value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(int value)
    {
        return init(static_cast<basetype>(value));
    }

    static always_inline type convert(floatx2_t value)
    {
        auto v = convert(static_cast<__m128>(value));
        return _mm256_shuffle_ps(v, v, _MM_SHUFFLE(1, 0, 1, 0));
    }

    static always_inline type convert(__m128 value)
    {
        return _mm256_set_m128(value, value);
    }
    static always_inline type convert(__m128i value)
    {
        return convert(intrin<float, 4>::convert(value));
    }
    static always_inline type convert(__m128d value)
    {
        return convert(intrin<float, 4>::convert(value));
    }

    static always_inline type convert(__m256i value)
    {
        return _mm256_cvtepi32_ps(value);
    }
    static always_inline type convert(__m256 value) { return value; }
    static always_inline type convert(__m256d value)
    {
        return convert(_mm256_cvtpd_ps(value));
    }
};

///////////////// DOUBLE x 4 ////////////////
template <> struct intrin<double, 4> {
    using basetype               = double;
    using type                   = __m256d;
    using maskbasetype           = basetype;
    using masktype               = type;
    static constexpr auto init   = _mm256_set1_pd;
    static constexpr auto load   = _mm256_load_pd;
    static constexpr auto loadu  = _mm256_loadu_pd;
    static constexpr auto store  = _mm256_store_pd;
    static constexpr auto storeu = _mm256_storeu_pd;
    static constexpr auto add    = _mm256_add_pd;
    static constexpr auto sub    = _mm256_sub_pd;
    static always_inline type vectorcall neg(type x)
    {
        return sub(init(basetype(0)), x);
    }
    static constexpr auto mul  = _mm256_mul_pd;
    static constexpr auto div  = _mm256_div_pd;
    static constexpr auto sqrt = _mm256_sqrt_pd;

    static always_inline masktype vectorcall signbit(type x)
    {
        auto shift = _mm256_srai_epi32((__m256i)x, 31);
        return (masktype)_mm256_shuffle_epi32(shift, _MM_SHUFFLE(3, 3, 1, 1));
    }

    static constexpr auto bitAnd = _mm256_and_pd;
    static constexpr auto bitOr  = _mm256_or_pd;
    static constexpr auto bitXor = _mm256_xor_pd;
    static always_inline type vectorcall bitNeg(type x)
    {
        return bitXor(init(basetype(0)), x);
    }
    static constexpr auto max = _mm256_max_pd;
    static constexpr auto min = _mm256_min_pd;

    static always_inline type vectorcall abs(type value)
    {
        return bitAnd(value, init(*reinterpret_cast<const basetype *>(
                                 &floatMask<basetype>::kNotSign)));
    }

    template <size_t Id, size_t K>
    static always_inline auto vectorcall getlane(type value)
    {
        if ((Id * K) % 2 == 1) {
            value = _mm256_shuffle_pd(value, value, _MM_SHUFFLE(0, 0, 1, 1));
        }

        auto r = _mm256_extractf128_pd(value, Id * K >= 2 ? 1 : 0);
        if constexpr (K == 1) {
            return _mm_cvtsd_f64(r);
        } else if constexpr (K == 2) {
            return r;
        } else
            return value;
    }

    template <int Shift> static type vectorcall shift(type value)
    {
        auto iv = _mm256_castpd_si256(value);
        auto vp = _mm256_permute2f128_si256(
            iv, iv, Shift > 0 ? (1 << 3) : 1 + (1 << 7));

        constexpr auto kShift = Shift > 0 ? Shift : -Shift;
        if constexpr (kShift == 2) return _mm256_castsi256_pd(vp);

        __m256i zero{};
        auto v1  = Shift > 0 ? iv : vp;
        auto v2  = Shift > 0 ? vp : iv;
        auto va1 = kShift < 2 ? v1 : v2;
        auto va2 = kShift < 2 ? v2 : zero;

        auto r = Shift == -3 ? _mm256_alignr_epi8(zero, vp, 8)
                             : _mm256_alignr_epi8(va1, va2, 8);

        return _mm256_castsi256_pd(r);
    }

    static always_inline type vectorcall flip1(type value)
    {
        return _mm256_shuffle_pd(value, value, _MM_SHUFFLE(0, 0, 1, 1));
    }
    static always_inline type vectorcall flip2(type value)
    {
        return _mm256_permute2f128_pd(value, value, _MM_SHUFFLE2(0, 3));
    }
    static always_inline basetype vectorcall sum(type value)
    {
        type flip = flip1(value);
        type sum  = add(value, flip);
        flip      = flip2(sum);
        sum       = add(sum, flip);
        return _mm256_cvtsd_f64(sum);
    }

    static always_inline __m128d vectorcall reduce2(type value)
    {
        auto flip = flip2(value);
        auto sum  = add(value, flip);
        return _mm256_extractf128_pd(sum, 0);
    }

    template <size_t K> static always_inline type vectorcall prefix(type value)
    {
        static_assert(K == 1 || K == 2);
        if (K == 1) {
            auto shift = _mm256_castsi256_pd(
                _mm256_slli_si256(_mm256_castpd_si256(value), 8));
            value = _mm256_add_pd(value, shift);
        }

        __m128d sum1block;
        switch (K) {
        default:
        case 1:
            sum1block = _mm_set1_pd(value[1]);
            break;
        case 2:
            sum1block = _mm256_castpd256_pd128(value);
            break;
        }
        return _mm256_add_pd(value, _mm256_setr_m128d(__m128d{}, sum1block));
    }

    static always_inline type push(type x, basetype v)
    {
        auto sh = shift<1>(x);
        auto lo = _mm_move_sd(_mm256_castpd256_pd128(sh), _mm_set_sd(v));
        return _mm256_insertf128_pd(sh, lo, 0);
    }
    static always_inline type push(type x, __m128d v)
    {
        auto sh = shift<2>(x);
        return _mm256_insertf128_pd(sh, v, 0);
    }

    static always_inline type vectorcall cmpeq(type x1, type x2)
    {
        return _mm256_cmp_pd(x1, x2, 0x00);
    }
    static always_inline type vectorcall cmpgt(type x1, type x2)
    {
        return _mm256_cmp_pd(x1, x2, 0x0E);
    }
    static always_inline type vectorcall cmpge(type x1, type x2)
    {
        return _mm256_cmp_pd(x1, x2, 0x0D);
    }
    static always_inline type vectorcall cmplt(type x1, type x2)
    {
        return _mm256_cmp_pd(x1, x2, 0x11);
    }
    static always_inline type vectorcall cmple(type x1, type x2)
    {
        return _mm256_cmp_pd(x1, x2, 0x12);
    }

    static always_inline type vectorcall blend(masktype m, type x2, type x1)
    {
        return _mm256_blendv_pd(x1, x2, m);
    }

    static constexpr auto any = _mm256_movemask_pd;
    static always_inline bool vectorcall all(masktype m)
    {
        return any(m) == 0x000f;
    }

    // convert
    static always_inline type convert(float value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(int value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(double value) { return init(value); }

    static always_inline type convert(floatx2_t value)
    {
        return convert(intrin<float, 4>::convert(value));
    }

    static always_inline type convert(__m128 value)
    {
        return _mm256_cvtps_pd(value);
    }
    static always_inline type convert(__m128i value)
    {
        return _mm256_cvtepi32_pd(value);
    }
    static always_inline type convert(__m128d value)
    {
        return _mm256_set_m128d(value, value);
    }

    static always_inline type convert(__m256i value)
    {
        return convert(intrin<int, 4>::convert(value));
    }
    static always_inline type convert(__m256 value)
    {
        return convert(intrin<float, 4>::convert(value));
    }
};
#endif

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp

#endif
