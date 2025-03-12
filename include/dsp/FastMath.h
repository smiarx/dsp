#pragma once

#include "Signal.h"
#include "Utils.h"
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

template <typename F, size_t N, size_t H = ilog2(N)>
Sample<F, N> hadamard(Sample<F, N> x)
{
    static constexpr auto L = 1 << H;

    // first we multiply values by power normalization
    if constexpr (N == L) {
        constexpr F powerNorm = std::pow(F(1) / F(2), F(H) / F(2));
        for (size_t i = 0; i < N; ++i) {
            x[i] *= powerNorm;
        }
    }

    if constexpr (H == 0) {
        return x;
    } else {
        static_assert(N % L == 0, "N must be divisible by L");
        constexpr auto L2 = L / 2;

        decltype(x) y;

#pragma omp simd
        for (size_t i = 0; i < N; i += L) {
            for (size_t j = 0; j < L2; ++j) {
                auto n    = i + j;
                y[n]      = x[n] + x[n + L2];
                y[n + L2] = x[n] - x[n + L2];
            }
        }
        return hadamard<F, N, H - 1>(y);
    }
}

} // namespace dsp
