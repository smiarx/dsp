#pragma once

#include "simd/multi.h"
#include <cmath>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

// constants
// NOLINTBEGIN (readability-identifier-naming)
template <typename T> struct constants {
    using F                       = baseType<T>;
    static constexpr auto pi      = F(3.14159265358979323846);
    static constexpr auto pi_2    = F(1.57079632679489661923);
    static constexpr auto pi_4    = F(0.78539816339744830962);
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

/* hermite interpolation */
template <typename F>
inline constexpr auto hermite(F ym1, F y0, F y1, F y2, F x)
{
    auto c0 = y0;
    auto c1 = F(.5) * (y1 - ym1);
    auto c2 = ym1 - F(2.5) * y0 + F(2.) * y1 - F(.5) * y2;
    auto c3 = F(.5) * (y2 - ym1) + F(1.5) * (y0 - y1);
    return c0 + (x * (c1 + x * (c2 + x * c3)));
}

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

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
