#pragma once

#if defined(__aarch64__)

#include "simd_default.h"
#include <arm_neon.h>

#define DSP_VEC_SIZE 4

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
    static constexpr auto div    = vdiv_f32;
    static constexpr auto sqrt   = vsqrt_f32;

    static always_inline auto vectorcall bitAnd(type x1, type x2)
    {
        return vand_u32(vreinterpret_u32_f32(x1), vreinterpret_u32_f32(x2));
    }
    static always_inline auto vectorcall bitOr(type x1, type x2)
    {
        return vorr_u32(vreinterpret_u32_f32(x1), vreinterpret_u32_f32(x2));
    }
    static always_inline auto vectorcall bitXor(type x1, type x2)
    {
        return veor_u32(vreinterpret_u32_f32(x1), vreinterpret_u32_f32(x2));
    }
    static always_inline auto vectorcall bitNeg(type x)
    {
        return bitXor(init(basetype(0)), x);
    }
    static constexpr auto max = vmax_f32;
    static constexpr auto min = vmin_f32;

    static constexpr auto abs = vabs_f32;
    static constexpr auto sum = vaddv_f32;

    static constexpr auto cmpeq = vceq_f32;
    static constexpr auto cmpgt = vcgt_f32;
    static constexpr auto blend = vbsl_f32;

    static always_inline auto vectorcall any(masktype x)
    {
        return vmaxv_u32(x) != 0;
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
    static always_inline type convert(float64x2_t value)
    {
        return vcvt_f32_f64(value);
    }
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
    static constexpr auto div    = vdivq_f32;
    static constexpr auto sqrt   = vsqrtq_f32;

    static always_inline auto vectorcall bitAnd(type x1, type x2)
    {
        return vandq_u32(vreinterpretq_u32_f32(x1), vreinterpretq_u32_f32(x2));
    }
    static always_inline auto vectorcall bitOr(type x1, type x2)
    {
        return vorrq_u32(vreinterpretq_u32_f32(x1), vreinterpretq_u32_f32(x2));
    }
    static always_inline auto vectorcall bitXor(type x1, type x2)
    {
        return veorq_u32(vreinterpretq_u32_f32(x1), vreinterpretq_u32_f32(x2));
    }
    static always_inline auto vectorcall bitNeg(type x)
    {
        return bitXor(init(basetype(0)), x);
    }
    static constexpr auto max = vmaxq_f32;
    static constexpr auto min = vminq_f32;

    static constexpr auto abs = vabsq_f32;
    static constexpr auto sum = vaddvq_f32;

    static constexpr auto cmpeq = vceqq_f32;
    static constexpr auto cmpgt = vcgtq_f32;
    static constexpr auto blend = vbslq_f32;

    static always_inline auto vectorcall any(masktype x)
    {
        return vmaxvq_u32(x) != 0;
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
    static always_inline type convert(int32x4_t value)
    {
        return vcvtq_f32_s32(value);
    }
    static always_inline type convert(float32x2_t value)
    {
        return vcombine_f32(value, value);
    }
    static always_inline type convert(float64x2_t value)
    {
        return convert(intrin<float, 2>::convert(value));
    }
};

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

    static always_inline auto vectorcall bitAnd(type x1, type x2)
    {
        return vandq_u64(vreinterpretq_u64_f64(x1), vreinterpretq_u64_f64(x2));
    }
    static always_inline auto vectorcall bitOr(type x1, type x2)
    {
        return vorrq_u64(vreinterpretq_u64_f64(x1), vreinterpretq_u64_f64(x2));
    }
    static always_inline auto vectorcall bitXor(type x1, type x2)
    {
        return veorq_u64(vreinterpretq_u64_f64(x1), vreinterpretq_u64_f64(x2));
    }
    static always_inline auto vectorcall bitNeg(type x)
    {
        return bitXor(init(basetype(0)), x);
    }
    static constexpr auto max = vmaxq_f64;
    static constexpr auto min = vminq_f64;

    static constexpr auto abs = vabsq_f64;
    static constexpr auto sum = vaddvq_f64;

    static constexpr auto cmpeq = vceqq_f64;
    static constexpr auto cmpgt = vcgtq_f64;
    static constexpr auto blend = vbslq_f64;

    static always_inline auto vectorcall any(masktype x)
    {
        return vmaxvq_f64(vreinterpretq_f64_u64(x)) != 0;
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

    static always_inline type convert(float32x2_t value)
    {
        return vcvt_f64_f32(value);
    }
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

    static always_inline auto vectorcall bitAnd(type x1, type x2)
    {
        return vandq_u32(vreinterpretq_u32_s32(x1), vreinterpretq_u32_s32(x2));
    }
    static always_inline auto vectorcall bitOr(type x1, type x2)
    {
        return vorrq_u32(vreinterpretq_u32_s32(x1), vreinterpretq_u32_s32(x2));
    }
    static always_inline auto vectorcall bitXor(type x1, type x2)
    {
        return veorq_u32(vreinterpretq_u32_s32(x1), vreinterpretq_u32_s32(x2));
    }
    static always_inline auto vectorcall bitNeg(type x)
    {
        return bitXor(init(basetype(0)), x);
    }
    static constexpr auto max = vmaxq_s32;
    static constexpr auto min = vminq_s32;

    static constexpr auto abs = vabsq_s32;
    static constexpr auto sum = vaddvq_s32;

    static constexpr auto cmpeq = vceqq_s32;
    static constexpr auto cmpgt = vcgtq_s32;
    static constexpr auto blend = vbslq_s32;

    static always_inline auto vectorcall any(masktype x)
    {
        return vmaxvq_u32(x) != 0;
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

    static always_inline type convert(float64x2_t value)
    {
        return convert(vcvt_f32_f64(value));
    }
};

} // namespace dsp
#endif
