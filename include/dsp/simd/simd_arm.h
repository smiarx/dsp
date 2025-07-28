#pragma once

#include "../cpu/defines.h"

#if DSP_NEON

#include "simd_default.h"
#include <arm_neon.h>

#define DSP_MAX_VEC_SIZE 16

#if DSP_AARCH64
#define DSP_SIMD_DOUBLE
#endif

namespace dsp
{

template <typename T, size_t N> struct intrin;

///////////////// FLOAT x 2 ////////////////
template <> struct intrin<float, 2> {
    using basetype               = float;
    using type                   = float32x2_t;
    using maskbasetype           = uint32_t;
    using masktype               = uint32x2_t;
    static constexpr auto init   = vdup_n_f32;
    static constexpr auto load   = vld1_f32;
    static constexpr auto loadu  = vld1_f32;
    static constexpr auto store  = vst1_f32;
    static constexpr auto storeu = vst1_f32;
    static constexpr auto add    = vadd_f32;
    static constexpr auto sub    = vsub_f32;
    static constexpr auto neg    = vneg_f32;
    static constexpr auto mul    = vmul_f32;
#ifdef DSP_AARCH64
    static constexpr auto div  = vdiv_f32;
    static constexpr auto sqrt = vsqrt_f32;
#else
    static always_inline type vectorcall div(type x, type y)
    {
        for (size_t i = 0; i < 2; ++i) x[i] /= y[i];
        return x;
    }
    static always_inline type vectorcall sqrt(type x)
    {
        for (size_t i = 0; i < 2; ++i) x[i] = std::sqrt(x[i]);
        return x;
    }
#endif

    static always_inline masktype vectorcall signbit(type x)
    {
        return (masktype)vshl_n_s32(vreinterpret_s32_f32(x), 31);
    }

    static constexpr auto bitAnd = vand_u32;
    static constexpr auto bitOr  = vorr_u32;
    static constexpr auto bitNeg = vmvn_u32;
    static always_inline masktype vectorcall bitXor(masktype x, masktype y)
    {
        return veor_u32(x, y);
    }
    static always_inline type vectorcall bitXor(type x, type y)
    {
        return vreinterpret_f32_u32(
            veor_u32(vreinterpret_u32_f32(x), vreinterpret_u32_f32(y)));
    }

    static constexpr auto max = vmax_f32;
    static constexpr auto min = vmin_f32;

    static constexpr auto abs = vabs_f32;
#ifdef DSP_AARCH64
    static constexpr auto sum = vaddv_f32;
#else
    static always_inline basetype vectorcall sum(type value)
    {
        return vget_lane_f32(vpadd_f32(value, value), 0);
    }
#endif

    static always_inline type vectorcall flip1(type value)
    {
        return vrev64_f32(value);
    }

    static constexpr auto cmpeq = vceq_f32;
    static constexpr auto cmpgt = vcgt_f32;
    static constexpr auto cmpge = vcge_f32;
    static constexpr auto cmplt = vclt_f32;
    static constexpr auto cmple = vcle_f32;
    static constexpr auto blend = vbsl_f32;

    static always_inline auto vectorcall any(masktype x)
    {
        return vget_lane_u64(vreinterpret_u64_u32(x), 0) != 0;
    }
    static always_inline auto vectorcall all(masktype x)
    {
        return vget_lane_u64(vreinterpret_u64_u32(x), 0) == 0xffffffffffffffff;
    }

    // convert
    static always_inline type convert(basetype value) { return init(value); }
    static always_inline type convert(int value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(double value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(int32x2_t value)
    {
        return vcvt_f32_s32(value);
    }
#ifdef DSP_AARCH64
    static always_inline type convert(float64x2_t value)
    {
        return vcvt_f32_f64(value);
    }
#endif
};

///////////////// FLOAT x 4 ////////////////
template <> struct intrin<float, 4> {
    using basetype               = float;
    using type                   = float32x4_t;
    using maskbasetype           = uint32_t;
    using masktype               = uint32x4_t;
    static constexpr auto init   = vdupq_n_f32;
    static constexpr auto load   = vld1q_f32;
    static constexpr auto loadu  = vld1q_f32;
    static constexpr auto store  = vst1q_f32;
    static constexpr auto storeu = vst1q_f32;
    static constexpr auto add    = vaddq_f32;
    static constexpr auto sub    = vsubq_f32;
    static constexpr auto neg    = vnegq_f32;
    static constexpr auto mul    = vmulq_f32;
#ifdef DSP_AARCH64
    static constexpr auto div  = vdivq_f32;
    static constexpr auto sqrt = vsqrtq_f32;
#else
    static always_inline type vectorcall div(type x, type y)
    {
        for (size_t i = 0; i < 4; ++i) x[i] /= y[i];
        return x;
    }
    static always_inline type vectorcall sqrt(type x)
    {
        for (size_t i = 0; i < 4; ++i) x[i] = std::sqrt(x[i]);
        return x;
    }
#endif

    static always_inline masktype vectorcall signbit(type x)
    {
        return (masktype)vshlq_n_s32(vreinterpretq_s32_f32(x), 31);
    }

    static constexpr auto bitAnd = vandq_u32;
    static constexpr auto bitOr  = vorrq_u32;
    static constexpr auto bitNeg = vmvnq_u32;

    static always_inline masktype vectorcall bitXor(masktype x, masktype y)
    {
        return veorq_u32(x, y);
    }
    static always_inline type vectorcall bitXor(type x, type y)
    {
        return vreinterpretq_f32_u32(
            veorq_u32(vreinterpretq_u32_f32(x), vreinterpretq_u32_f32(y)));
    }

    static constexpr auto max = vmaxq_f32;
    static constexpr auto min = vminq_f32;

    static constexpr auto abs = vabsq_f32;
#ifdef DSP_AARCH64
    static constexpr auto sum = vaddvq_f32;
#else
    static always_inline basetype vectorcall sum(type value)
    {
        auto sum = vpadd_f32(vget_low_f32(value), vget_high_f32(value));
        return vget_lane_f32(vpadd_f32(sum, sum), 0);
    }
#endif

    static always_inline type vectorcall flip1(type value)
    {
        return vrev64q_f32(value);
    }
    static always_inline type vectorcall flip2(type value)
    {
#ifdef DSP_AARCH64
        auto valued = vreinterpretq_f64_f32(value);
        auto shufd  = vextq_f64(valued, valued, 1);
        return vreinterpretq_f32_f64(shufd);
#else
        return vcombine_f32(vget_high_f32(value), vget_low_f32(value));
#endif
    }
    static always_inline float32x2_t vectorcall reduce2(type x)
    {
        return vget_low_f32(add(x, flip2(x)));
    }

    static constexpr auto cmpeq = vceqq_f32;
    static constexpr auto cmpgt = vcgtq_f32;
    static constexpr auto cmpge = vcgeq_f32;
    static constexpr auto cmplt = vcltq_f32;
    static constexpr auto cmple = vcleq_f32;
    static constexpr auto blend = vbslq_f32;

    static always_inline auto vectorcall any(masktype x)
    {
        return vget_lane_u64(vreinterpret_u64_u16(vqmovn_u32(x)), 0) != 0;
    }
    static always_inline auto vectorcall all(masktype x)
    {
        return vget_lane_u64(vreinterpret_u64_u16(vqmovn_u32(x)), 0) ==
               0xffffffffffffffff;
    }

    // convert
    static always_inline type convert(basetype value) { return init(value); }
    static always_inline type convert(int value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(double value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(int32x2_t value)
    {
        return convert(vcvt_f32_s32(value));
    }
    static always_inline type convert(int32x4_t value)
    {
        return vcvtq_f32_s32(value);
    }
    static always_inline type convert(float32x2_t value)
    {
        return vcombine_f32(value, value);
    }
#ifdef DSP_AARCH64
    static always_inline type convert(float64x2_t value)
    {
        return convert(intrin<float, 2>::convert(value));
    }
#endif
};

#ifdef DSP_AARCH64
///////////////// DOUBLE x 2 ////////////////
template <> struct intrin<double, 2> {
    using basetype               = double;
    using type                   = float64x2_t;
    using maskbasetype           = uint64_t;
    using masktype               = uint64x2_t;
    static constexpr auto init   = vdupq_n_f64;
    static constexpr auto load   = vld1q_f64;
    static constexpr auto loadu  = vld1q_f64;
    static constexpr auto store  = vst1q_f64;
    static constexpr auto storeu = vst1q_f64;
    static constexpr auto add    = vaddq_f64;
    static constexpr auto sub    = vsubq_f64;
    static constexpr auto neg    = vnegq_f64;
    static constexpr auto mul    = vmulq_f64;
    static constexpr auto div    = vdivq_f64;
    static constexpr auto sqrt   = vsqrtq_f64;

    static always_inline masktype vectorcall signbit(type x)
    {
        return (masktype)vshlq_n_s64(vreinterpretq_s64_f64(x), 31);
    }

    static constexpr auto bitAnd = vandq_u64;
    static constexpr auto bitOr  = vorrq_u64;
    static always_inline masktype vectorcall bitXor(masktype x, masktype y)
    {
        return veorq_u64(x, y);
    }
    static always_inline type vectorcall bitXor(type x, type y)
    {
        return vreinterpretq_f64_u64(
            veorq_u64(vreinterpretq_u64_f64(x), vreinterpretq_u64_f64(y)));
    }
    static always_inline auto vectorcall bitNeg(masktype x)
    {
        bitXor(vdupq_n_u64(0), x);
    }

    static constexpr auto max = vmaxq_f64;
    static constexpr auto min = vminq_f64;

    static constexpr auto abs = vabsq_f64;

    static constexpr auto sum = vaddvq_f64;
    static always_inline type vectorcall flip1(type value)
    {
        return vextq_f64(value, value, 1);
    }

    static constexpr auto cmpeq = vceqq_f64;
    static constexpr auto cmpgt = vcgtq_f64;
    static constexpr auto cmpge = vcgeq_f64;
    static constexpr auto cmplt = vcltq_f64;
    static constexpr auto cmple = vcleq_f64;
    static constexpr auto blend = vbslq_f64;

    static always_inline auto vectorcall any(masktype x)
    {
        return vget_lane_u64(vreinterpret_u64_u32(vqmovn_u64(x)), 0) != 0;
    }
    static always_inline auto vectorcall all(masktype x)
    {
        return vget_lane_u64(vreinterpret_u64_u32(vqmovn_u64(x)), 0) ==
               0xffffffffffffffff;
    }

    // convert
    static always_inline type convert(basetype value) { return init(value); }
    static always_inline type convert(int value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(float value)
    {
        return init(static_cast<basetype>(value));
    }

    static always_inline type convert(int32x2_t value)
    {
        return vcvt_f64_f32(vcvt_f32_s32(value));
    }

    static always_inline type convert(float32x2_t value)
    {
        return vcvt_f64_f32(value);
    }
};
#endif

///////////////// INT x 2 ////////////////
template <> struct intrin<int32_t, 2> {
    using basetype               = int32_t;
    using type                   = int32x2_t;
    using maskbasetype           = uint32_t;
    using masktype               = uint32x2_t;
    static constexpr auto init   = vdup_n_s32;
    static constexpr auto load   = vld1_s32;
    static constexpr auto loadu  = vld1_s32;
    static constexpr auto store  = vst1_s32;
    static constexpr auto storeu = vst1_s32;
    static constexpr auto add    = vadd_s32;
    static constexpr auto sub    = vsub_s32;
    static constexpr auto neg    = vneg_s32;
    static constexpr auto mul    = vmul_s32;

    static always_inline auto vectorcall bitAnd(masktype x1, masktype x2)
    {
        return vand_u32(x1, x2);
    }
    static always_inline auto vectorcall bitOr(masktype x1, masktype x2)
    {
        return vorr_u32(x1, x2);
    }
    static always_inline auto vectorcall bitXor(masktype x1, masktype x2)
    {
        return veor_u32(x1, x2);
    }
    static always_inline auto vectorcall bitNeg(masktype x)
    {
        return vmvn_u32(x);
    }

    static always_inline auto vectorcall bitAnd(type x1, type x2)
    {
        return vand_s32(x1, x2);
    }
    static always_inline auto vectorcall bitOr(type x1, type x2)
    {
        return vorr_s32(x1, x2);
    }
    static always_inline auto vectorcall bitXor(type x1, type x2)
    {
        return veor_s32(x1, x2);
    }
    static always_inline auto vectorcall bitNeg(type x) { return vmvn_s32(x); }

    static constexpr auto bitShiftLeft  = vshl_n_s32;
    static constexpr auto bitShiftRight = vshr_n_s32;

    static constexpr auto max = vmax_s32;
    static constexpr auto min = vmin_s32;

    static constexpr auto abs = vabs_s32;
#ifdef DSP_AARCH64
    static constexpr auto sum = vaddv_s32;
#else
    static always_inline basetype vectorcall sum(type value)
    {
        return vget_lane_s32(vpadd_s32(value, value), 0);
    }
#endif

    static always_inline type vectorcall flip1(type value)
    {
        return vrev64_s32(value);
    }

    static constexpr auto cmpeq = vceq_s32;
    static constexpr auto cmpgt = vcgt_s32;
    static constexpr auto cmpge = vcge_s32;
    static constexpr auto cmplt = vclt_s32;
    static constexpr auto cmple = vcle_s32;
    static constexpr auto blend = vbsl_s32;

    static always_inline auto vectorcall any(masktype x)
    {
#ifdef DSP_AARCH64
        return vmaxv_u32(x) != 0;
#else
        return vget_lane_u64(vreinterpret_u64_u32(x), 0) != 0;
#endif
    }
    static always_inline auto vectorcall all(masktype x)
    {
#ifdef DSP_AARCH64
        return vminv_u32(x) != 0;
#else
        return vget_lane_u64(vreinterpret_u64_u32(x), 0) == 0xffffffffffffffff;
#endif
    }

    // convert
    static always_inline type convert(basetype value) { return init(value); }
    static always_inline type convert(double value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(float value)
    {
        return init(static_cast<basetype>(value));
    }

    static always_inline type convert(float32x2_t value)
    {
        return vcvt_s32_f32(value);
    }

#ifdef DSP_AARCH64
    static always_inline type convert(float64x2_t value)
    {
        return convert(vcvt_f32_f64(value));
    }
#endif
};

///////////////// INT x 4 ////////////////
template <> struct intrin<int32_t, 4> {
    using basetype               = int32_t;
    using type                   = int32x4_t;
    using maskbasetype           = uint32_t;
    using masktype               = uint32x4_t;
    static constexpr auto init   = vdupq_n_s32;
    static constexpr auto load   = vld1q_s32;
    static constexpr auto loadu  = vld1q_s32;
    static constexpr auto store  = vst1q_s32;
    static constexpr auto storeu = vst1q_s32;
    static constexpr auto add    = vaddq_s32;
    static constexpr auto sub    = vsubq_s32;
    static constexpr auto neg    = vnegq_s32;
    static constexpr auto mul    = vmulq_s32;

    static always_inline auto vectorcall bitAnd(masktype x1, masktype x2)
    {
        return vandq_u32(x1, x2);
    }
    static always_inline auto vectorcall bitOr(masktype x1, masktype x2)
    {
        return vorrq_u32(x1, x2);
    }
    static always_inline auto vectorcall bitXor(masktype x1, masktype x2)
    {
        return veorq_u32(x1, x2);
    }
    static always_inline auto vectorcall bitNeg(masktype x)
    {
        return vmvnq_u32(x);
    }

    static always_inline auto vectorcall bitAnd(type x1, type x2)
    {
        return vandq_s32(x1, x2);
    }
    static always_inline auto vectorcall bitOr(type x1, type x2)
    {
        return vorrq_s32(x1, x2);
    }
    static always_inline auto vectorcall bitXor(type x1, type x2)
    {
        return veorq_s32(x1, x2);
    }
    static always_inline auto vectorcall bitNeg(type x) { return vmvnq_s32(x); }

    static constexpr auto bitShiftLeft  = vshlq_n_s32;
    static constexpr auto bitShiftRight = vshrq_n_s32;

    static constexpr auto max = vmaxq_s32;
    static constexpr auto min = vminq_s32;

    static constexpr auto abs = vabsq_s32;

#ifdef DSP_AARCH64
    static constexpr auto sum = vaddvq_s32;
#else
    static always_inline basetype vectorcall sum(type value)
    {
        auto sum = vpadd_s32(vget_low_s32(value), vget_high_s32(value));
        return vget_lane_s32(vpadd_s32(sum, sum), 0);
    }
#endif

    static always_inline type vectorcall flip1(type value)
    {
        return vrev64q_s32(value);
    }
    static always_inline type vectorcall flip2(type value)
    {
        auto valued = vreinterpretq_s64_s32(value);
        auto shufd  = vextq_s64(valued, valued, 1);
        return vreinterpretq_s32_s64(shufd);
    }

    static always_inline auto vectorcall reduce2(type x)
    {
        return vget_low_s32(add(x, flip2(x)));
    }

    static constexpr auto cmpeq = vceqq_s32;
    static constexpr auto cmpgt = vcgtq_s32;
    static constexpr auto cmpge = vcgeq_s32;
    static constexpr auto cmplt = vcltq_s32;
    static constexpr auto cmple = vcleq_s32;
    static constexpr auto blend = vbslq_s32;

    static always_inline auto vectorcall any(masktype x)
    {
        return vget_lane_u64(vreinterpret_u64_u16(vqmovn_u32(x)), 0) != 0;
    }
    static always_inline auto vectorcall all(masktype x)
    {
        return vget_lane_u64(vreinterpret_u64_u16(vqmovn_u32(x)), 0) ==
               0xffffffffffffffff;
    }

    // convert
    static always_inline type convert(basetype value) { return init(value); }
    static always_inline type convert(double value)
    {
        return init(static_cast<basetype>(value));
    }
    static always_inline type convert(float value)
    {
        return init(static_cast<basetype>(value));
    }

    static always_inline type convert(float32x4_t value)
    {
        return vcvtq_s32_f32(value);
    }

    static always_inline type convert(int32x2_t value)
    {
        return vcombine_s32(value, value);
    }

    static always_inline type convert(float32x2_t value)
    {
        return convert(vcvt_s32_f32(value));
    }

#ifdef DSP_AARCH64
    static always_inline type convert(float64x2_t value)
    {
        return convert(vcvt_f32_f64(value));
    }
#endif
};

} // namespace dsp
#endif
