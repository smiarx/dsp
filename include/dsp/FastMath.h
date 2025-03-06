#pragma once

#include <cmath>

namespace dsp
{

// https://varietyofsound.wordpress.com/2011/02/14/efficient-tanh-computation-using-lamberts-continued-fraction/
template <typename F> static inline constexpr F tanh(F x)
{
    F xsq    = x * x;
    F num    = x * (F(135135) + xsq * (F(17325) + xsq * (F(378) + xsq)));
    F den    = F(135135.) + xsq * (F(62370) + xsq * (F(3150) + F(28) * xsq));
    F result = num / den;
    result   = result > F(1) ? F(1) : result;
    result   = result < -F(1) ? -F(1) : result;
    return result;
}

template <typename Float> constexpr auto sinc(Float x)
{
    auto xpi = x * Float(M_PI);
    return fabs(x) < 0.0001f ? 1.f : sinf(xpi) / (xpi);
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

} // namespace dsp
