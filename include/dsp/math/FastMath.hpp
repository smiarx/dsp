#pragma once

#include <cassert>
#include <cmath>
#include <dsp/simd/multi.hpp>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

// constants
// NOLINTBEGIN (readability-identifier-naming)
template <typename T> struct constants {
    using F                       = baseType<T>;
    static constexpr auto twopi   = F(6.28318530717958647693);
    static constexpr auto pi      = F(3.14159265358979323846);
    static constexpr auto pi_2    = F(1.57079632679489661923);
    static constexpr auto pi_4    = F(0.78539816339744830962);
    static constexpr auto invPi   = F(0.31830988618379067154);
    static constexpr auto sqrt1_2 = F(0.70710678118654752440);
};

// NOLINTEND (readability-identifier-naming)

// https://varietyofsound.wordpress.com/2011/02/14/efficient-tanh-computation-using-lamberts-continued-fraction/
// https://math.stackexchange.com/questions/107292/rapid-approximation-of-tanhx
template <typename T> static inline constexpr T fasttanh(T x)
{
    /* if tanh(z/2) = num/den, then
     * tan(z) = 2⋅tanh(z/2)/(1+tanh²(z/2))
     *        = 2⋅num/den / (1+num²/den²)
     *        = 2⋅num⋅den / (num²+den²)
     */
    x     = x * T(0.5);
    T xsq = x * x;
    T num, den;
    if constexpr (std::is_same_v<baseType<T>, double>) {
        num = T(21) * x * (T(495) + xsq * (T(60) + xsq));
        den = T(10395) + xsq * (T(4725) + xsq * (T(210) + xsq));
    } else {
        num = x * (T(945) + xsq * (T(105) + xsq));
        den = T(15) * (T(63) + xsq * (T(28) + xsq));
    }
    T tanh = T(2) * num * den / (num * num + den * den);
    return tanh;
}

template <typename F> constexpr auto sinc(F x)
{
    auto xpi = x * constants<F>::pi;
    return fabs(x) < F(0.0001) ? F(1) : std::sin(xpi) / (xpi);
}

/* 0th order bessel value */
template <typename F> static constexpr auto zerothOrderBessel(F x)
{
    const auto eps = F(0.000001);

    //  initialize the series term for m=0 and the result
    F bessel = 0;
    F term   = 1;
    F m      = 0;

    //  accumulate terms as long as they are significant
    while (term > eps * bessel) {
        bessel += term;

        //  update the term
        m += F(1);
        term *= (x * x) / (4 * m * m);
    }

    return bessel;
}

///////////////////////////////////////////////////////////////

template <typename T> T fastexp2(T x)
{
    using bT                 = baseType<T>;
    constexpr bool kIsDouble = std::is_same_v<double, bT>;
    constexpr auto kFPMinExpo =
        bT(kIsDouble ? -1022 : -126); // exponent of minimum binary32 normal
    constexpr auto kFPExpoBias =
        bT(kIsDouble ? 1023 : 127); // binary32 exponent bias
    constexpr int kFPMantBits =
        kIsDouble ? 52 : 23; // number of stored mantissa (significand) bits
    constexpr auto kFPMult = bT((1ll << kFPMantBits));

    x = blend(x < kFPMinExpo, T(kFPMinExpo), x); // clamp
                                                 //
    /* 2**p = 2**(w+z), with w an integer and z in [0, 1) */
    auto w = floor(x); // integral part
    auto z = x - w;    // fractional part
    /* approximate 2**z-1 for z in [0, 1) */

    T approx;
    if constexpr (kIsDouble) {
        approx = T(-0x1.6e75d58p+2) +
                 T(0x1.bba7414p+4) / (T(0x1.35eccbap+2) - z) -
                 T(0x1.f5e53c2p-2) * z;
    } else {
        approx = T(-0x1.6e7592p+2f) +
                 T(0x1.bba764p+4f) / (T(0x1.35ed00p+2f) - z) -
                 T(0x1.f5e546p-2f) * z;
    }

    /* assemble the exponent and mantissa components into final result */
    auto res = kFPMult * (w + kFPExpoBias + approx);

    if constexpr (kIsDouble) {
        // no double simd right now...
        auto ires = static_cast<int64_t>(res);
        return std::bit_cast<T>(ires);
    } else {
        auto ires = toInt(res);
        return std::bit_cast<T>(ires);
    }
}

template <typename T> T fastlog2(T x)
{
    assert(all(x > T(0)));

    using bT                 = baseType<T>;
    constexpr bool kIsDouble = std::is_same_v<double, bT>;
    using iT = std::conditional_t<kIsDouble, int64_t, intType<T>>;
    constexpr auto kFPExpoBias = kIsDouble ? 1023l : 127l; // exponent bias
    constexpr int kFPMantBits =
        kIsDouble ? 52 : 23; // number of stored mantissa (significand) bits
    constexpr auto kFPMantMask = (1l << kFPMantBits) - 1;

    auto ix = std::bit_cast<iT>(x);

    // no need to mask sign bit, negative values should crash
    auto exp = (ix >> kFPMantBits) - iT(kFPExpoBias);

    // change exponent so that value is in [1,2[
    ix = iT(kFPExpoBias << kFPMantBits) | (iT(kFPMantMask) & ix);
    x  = std::bit_cast<T>(ix);

    auto approx = T(0x1.c8dfbep+0f) -
                  T(0x1.bd56d4p+1f) / (x + T(0x1.6bc5c6p-1f)) +
                  x * T(0x1.ff08a2p-3f);

    if constexpr (kIsDouble) {
        return static_cast<double>(exp) + approx;
    } else {
        return toFloat<bT>(exp) + approx;
    }
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
