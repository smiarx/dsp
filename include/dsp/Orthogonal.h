#pragma once

#include "FastMath.h"
#include "LinAlg.h"
#include "Utils.h"
#include <cmath>

namespace dsp
{

/* orthogonal transforms */

template <typename F, size_t N> Sample<F, N> householder(Sample<F, N> x)
{
    // householder transform

    F sum{};
    for (size_t i = 0; i < N; ++i) {
        sum += x[i];
    }
    sum /= F(N) / F(2);

    Sample<F, N> y;
    for (size_t i = 0; i < N; ++i) {
        y[i] = x[i] - sum;
    }

    return y;
}

template <typename F, size_t N, size_t H = ilog2(N)>
Sample<F, N> hadamardBase(Sample<F, N> x)
{
    /* compute hadamard transform of sample (without power normalization) */
    static constexpr auto kL = 1 << H;

    if constexpr (H == 0) {
        return x;
    } else {
        static_assert(N % kL == 0, "N must be divisible by L");
        constexpr auto kL2 = kL / 2;

        decltype(x) y;

#pragma omp simd
        for (size_t i = 0; i < N; i += kL) {
            for (size_t j = 0; j < kL2; ++j) {
                auto n     = i + j;
                y[n]       = x[n] + x[n + kL2];
                y[n + kL2] = x[n] - x[n + kL2];
            }
        }
        return hadamardBase<F, N, H - 1>(y);
    }
}

template <typename F, size_t N, size_t H = ilog2(N)>
Sample<F, N> hadamard(Sample<F, N> x)
{
    /* hadamard transform with power normalization */

    // first we multiply values by power normalization
    const F kPowerNorm = std::pow(F(1) / F(2), F(H) / F(2));
    for (size_t i = 0; i < N; ++i) {
        x[i] *= kPowerNorm;
    }

    return hadamardBase<F, N, H>(x);
}

template <size_t N>
static inline linalg::fMatrix<N> hadamardInterpolMatrix(float t)
{
    /*
     *  We want to interpolate between identity matrix and hadamard matrix
     *  ie a matrix H_(t) such as H_(0) = Id, H_(1) = H and H_(t) is orthogonal
     *  for every t in [0,1]
     *
     *  https://dafx2020.mdw.ac.at/proceedings/papers/DAFx2020_paper_53.pdf
     *  A solution is
     *
     *  H_(t) = exp(t⋅log(H))
     *
     *  We can diagonalise H = UDU^T, then
     *  H_(t) = Uexp(t⋅log(D))U^T
     *
     *      | 1  1  1  1 |
     *      | 1 -1  1 -1 |
     * H=1/2| 1  1 -1 -1 |
     *      | 1 -1 -1  1 |
     *
     * 4 eigen vectors
     *     v1 = |3  1  1  1| / sqrt(12) eigen value 1
     *     v2 = |0  1  1 -2| / sqrt(6)  eigen value 1
     *     v3 = |1 -1 -1 -1| / 2        eigen value -1
     *     v4 = |0  1 -1  0| / sqrt(2)  eigen value -1
     *
     *     |1  0  0  0|
     * D = |0  1  0  0|  and U = |v1 v2 v3 v4|
     *     |0  0 -1  0|
     *     |0  0  0 -1|
     *
     *
     *          |0  0  0   0|                  |1  0    0          0      |
     * log(D) = |0  0  0   0|  exp(t⋅log(D)) = |0  1    0          0      |
     *          |0  0  0 -pi|                  |0  0 cos(t⋅pi) -sin(t⋅pi) |
     *          |0  0 pi   0|                  |0  0 sin(t⋅pi)  cos(t⋅pi) |
     *
     * see
     * https://www.m-hikari.com/ija/ija-password-2008/ija-password1-4-2008/morsyIJA1-4-2008.pdf
     *
     * with
     *                        | 3/sqrt(12)     0      |
     * A = |v1 v2||v1 v2|^T = | 1/sqrt(12)  1/sqrt(6) | | 3/sqrt(12) 1/sqrt(12)
     * 1/sqrt(12) 1/sqrt(12) | | 1/sqrt(12)  1/sqrt(6) | |     0      1/sqrt(6)
     * 1/sqrt(6)  -2/sqrt(6) | | 1/sqrt(12) -2/sqrt(6) |
     *
     *                        | 3/4  1/4  1/4  1/4 |
     *                      = | 1/4  1/4  1/4 -1/4 |
     *                        | 1/4  1/4  1/4 -1/4 |
     *                        | 1/4 -1/4 -1/4  3/4 |
     *
     *
     *                        |  1/2  0        |
     * B = |v3 v4||v3 v4|^T = | -1/2  1/sqrt(2)|| 1/2   -1/2       -1/2    -1/2|
     *                        | -1/2 -1/sqrt(2)|| 0   1/sqrt(2) -1/sqrt(2)    0|
     *                        | -1/2  0        |
     *
     *                        |  1/4 -1/4 -1/4 -1/4 |
     *                        | -1/4  3/4 -1/4  1/4 |
     *                        | -1/4 -1/4  3/4  1/4 |
     *                        | -1/4  1/4  1/4  1/4 |
     *
     *                         |  1/2  0        |
     * C = |v3 v4||-v4 v3|^T = | -1/2  1/sqrt(2)|| 0   -1/sqrt(2)  1/sqrt(2) 0|
     *                         | -1/2 -1/sqrt(2)|| 1/2   -1/2       -1/2 -1/2|
     *                         | -1/2  0        |
     *
     *                         |     0      -1/sqrt(8)   1/sqrt(8)   0        |
     *                         | 1/sqrt(8)     0        -1/sqrt(2)  -1/sqrt(8)|
     *                         |-1/sqrt(8)   1/sqrt(2)     0         1/sqrt(8)|
     *                         |     0       1/sqrt(8)  -1/sqrt(8)    0       |
     *
     *
     * H_(t) = A + cos(t⋅pi)⋅B + sin(t⋅pi)⋅C
     *
     */

    static_assert(N == 4, "Implemented only for N=4");
    linalg::fMatrix<N> matA = {{{
        {0.75f, 0.25f, 0.25f, 0.25f},
        {0.25f, 0.25f, 0.25f, -0.25f},
        {0.25f, 0.25f, 0.25f, -0.25f},
        {0.25f, -0.25f, -0.25f, 0.75f},
    }}};
    linalg::fMatrix<N> matB = {{{
        {0.25f, -0.25f, -0.25f, -0.25f},
        {-0.25f, 0.75f, -0.25f, 0.25f},
        {-0.25f, -0.25f, 0.75f, 0.25f},
        {-0.25f, 0.25f, 0.25f, 0.25f},
    }}};

    constexpr auto kV12     = constants<float>::sqrt1_2;
    constexpr auto kV18     = 0.5f * kV12;
    linalg::fMatrix<N> matC = {{{
        {0.f, kV18, -kV18, 0.f},
        {-kV18, 0.f, kV12, kV18},
        {kV18, -kV12, 0.f, -kV18},
        {0.f, -kV18, kV18, 0.f},
    }}};

    linalg::fMatrix<N> matH{};

    auto c = cosf(dsp::constants<float>::pi * t);
    auto s = sinf(dsp::constants<float>::pi * t);
#pragma omp simd
    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < N; ++i) {
            matH[j][i] = matA[j][i] + c * matB[j][i] + s * matC[j][i];
        }
    }

    return matH;
}

} // namespace dsp
