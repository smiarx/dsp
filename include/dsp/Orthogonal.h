#pragma once

#include "FastMath.h"
#include "LinAlg.h"
#include "MultiVal.h"
#include "Utils.h"
#include <cmath>
#include <iostream>

namespace dsp
{

/* orthogonal transforms */

template <typename T> auto householder(T x)
{
    // householder transform
    auto s = sum(x);
    s /= kTypeWidth<T> / 2;

    return x - s;
}

template <typename T, size_t H = kTypeWidth<T>> auto hadamardBase(T x)
{
    // use matrix decomposition
    // | 1  1       || 1     1    |  | 1  1  1  1 |
    // | 1 -1       ||    1     1 |  | 1 -1  1 -1 |
    // |       1  1 || 1    -1    | =| 1  1 -1 -1 |
    // |       1 -1 ||    1    -1 |  | 1 -1 -1  1 |
    using namespace loadfuncs;

    /* compute hadamard transform of sample (without power normalization) */
    static constexpr auto kN = kTypeWidth<T>;

    if constexpr (H == 1) {
        return x;
    } else {
        static_assert(kN % H == 0, "N must be divisible by H");
        constexpr auto kH2 = H / 2;

        // flip
        auto xFlip = flip<kH2>(x);

        // reverse signs
        constexpr auto kSignMask = floatMask<baseType<T>>::kSign;
        MultiVal<baseType<T>, kN> mask;
        for (size_t i = 0; i < kN; ++i) {
            mask[i] = i & kH2 ? kSignMask.f : 0;
        }
        auto xNeg = load(mask) ^ load(x);

        auto y = xNeg + xFlip;
        return hadamardBase<T, kH2>(y);
    }
}
template <typename T, size_t H = kTypeWidth<T>> auto hadamard(T x)
{
    using namespace loadfuncs;
    /* hadamard transform with power normalization */

    // first we multiply values by power normalization
    if constexpr (H > 1) {
        if constexpr (H == 2) x = load(x) * constants<baseType<T>>::sqrt1_2;
        else if constexpr (H == 4)
            x = load(x) * 0.5;
        else if constexpr (H == 8)
            x = load(x) * 0.3535533905932738;
    }

    return hadamardBase<T, H>(x);
}

// inline float hadamard(float x) { return x;}
// inline double hadamard(double x) { return x;}
//
//// hadamardx2
//// | 1  1 |
//// | 1 -1 |
////
// #if defined(__x86_64__)
// inline simd<float,2> hadamard(floatx2_t x)
//{
//     constexpr auto kMask = floatMask<float>::kSign;
//
//     // x * sqrt(1/2)
//     x = _mm_mul_ps(x, _mm_set1_ps(constants<float>::sqrt1_2));
//
//     auto shuf = _mm_shuffle_ps(x,x,_MM_SHUFFLE(2,3,0,1));
//     auto m = _mm_set_ps(kMask.f,0,kMask.f,0);
//     x = _mm_xor_ps(x,m);
//
//     x = _mm_add_ps(x,shuf);
//
//     return floatx2_t(x);
// }
// inline simd<double,2> hadamard(__m128d x)
//{
//     constexpr auto kMask = floatMask<double>::kSign;
//
//     // x * sqrt(1/2)
//     x = _mm_mul_pd(x, _mm_set1_pd(constants<double>::sqrt1_2));
//
//     auto shuf = _mm_shuffle_pd(x,x,_MM_SHUFFLE2(0,1));
//     auto m = _mm_set_pd(kMask.f,0);
//     x = _mm_xor_pd(x,m);
//
//     return _mm_add_pd(x,shuf);
// }
// #elif defined(__aarch64__)
// inline simd<float,2> hadamard(float32x2_t x)
//{
//     constexpr auto kMask = floatMask<float>::kSign;
//
//     x = vmul_f32(x, vdup_n_f32(constants<float>::sqrt1_2));
//
//     auto shuf = vrev64_f32(x);
//     uint32_t m[] = {0,kMask.i};
//     auto vm = vld1_u32(m);
//     x = vreinterpret_f32_u32(veor_u32(vreinterpret_u32_f32(x),vm));
//
//     x = vadd_f32(x,shuf);
//     return x;
// }
// inline simd<double,2> hadamard(float64x2_t x)
//{
//     constexpr auto kMask = floatMask<double>::kSign;
//
//     x = vmulq_f64(x, vdupq_n_f64(constants<double>::sqrt1_2));
//
//     auto shuf = vextq_f64(x,x,1);
//     uint64_t m[] = {0,kMask.i};
//     auto vm = vld1q_u64(m);
//     x = vreinterpretq_f64_u64(veorq_u64(vreinterpretq_u64_f64(x),vm));
//
//     x = vaddq_f64(x,shuf);
//     return x;
// }
// #endif
//
//// hadamardx4
//// | 1     1    || 1  1       |  | 1  1  1  1 |
//// |    1     1 || 1 -1       |  | 1 -1  1 -1 |
//// | 1    -1    ||       1  1 | =| 1  1 -1 -1 |
//// |    1    -1 ||       1 -1 |  | 1 -1 -1  1 |
////
//
// #if defined(__x86_64__)
// inline simd<float,4> hadamard(__m128 x)
//{
//    constexpr auto kMask = floatMask<float>::kSign;
//
//    // x * sqrt(1/4)
//    x = _mm_mul_ps(x, _mm_set1_ps(0.5f));
//
//    auto shuf = _mm_shuffle_ps(x,x,_MM_SHUFFLE(2,3,0,1));
//    auto m = _mm_set_ps(kMask.f,0,kMask.f,0);
//    x = _mm_xor_ps(x,m);
//
//    x = _mm_add_ps(x,shuf);
//
//    shuf = _mm_shuffle_ps(x,x,_MM_SHUFFLE(1,0,3,2));
//    m = _mm_set_ps(kMask.f,kMask.f,0,0);
//    x = _mm_xor_ps(x,m);
//
//    x = _mm_add_ps(x,shuf);
//
//    return x;
//}
// #if defined(__AVX__)
// inline simd<double,4> hadamard(__m256d x)
//{
//    constexpr auto kMask = floatMask<double>::kSign;
//
//    // x * sqrt(1/4)
//    x = _mm256_mul_pd(x, _mm256_set1_pd(0.5));
//
//    auto shuf = _mm256_shuffle_pd(x,x, (0 << 3) + (1<<2) + (0 << 1) + 1);
//    auto m = _mm256_set_pd(kMask.f,0,kMask.f,0);
//    x = _mm256_xor_pd(x,m);
//
//    x = _mm256_add_pd(x,shuf);
//
//    shuf = _mm256_permute2f128_pd(x,x,1);
//    m = _mm256_set_pd(kMask.f,kMask.f,0,0);
//    x = _mm256_xor_pd(x,m);
//
//    x = _mm256_add_pd(x,shuf);
//
//    return x;
//}
// #endif
// #elif defined(__aarch64__)
// inline simd<float,4> hadamard(float32x4_t x)
//{
//    constexpr auto kMask = floatMask<float>::kSign;
//
//    // x * sqrt(1/4)
//    x = vmulq_f32(x, vdupq_n_f32(0.5f));
//    std::cout << "x\t" << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << x[3]
//    << "\n";
//
//    auto shuf = vrev64q_f32(x);
//    std::cout << "shuf\t" << shuf[0] << "\t" << shuf[1] << "\t" << shuf[2] <<
//    "\t" << shuf[3] << "\n";
//    {
//    uint32_t m[] = {0,kMask.i,0,kMask.i};
//    auto vm = vld1q_u32(m);
//    x = vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(x),vm));
//    std::cout << "x\t" << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << x[3]
//    << "\n";
//    }
//
//    x = vaddq_f32(x,shuf);
//    std::cout << "x\t" << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << x[3]
//    << "\n";
//
//    {
//    auto xd = vreinterpretq_f64_f32(x);
//    auto shufd = vextq_f64(xd,xd,1);
//    shuf = vreinterpretq_f32_f64(shufd);
//    }
//    std::cout << "shuf\t" << shuf[0] << "\t" << shuf[1] << "\t" << shuf[2] <<
//    "\t" << shuf[3] << "\n";
//    {
//    uint32_t m[] = {0,0,kMask.i,kMask.i};
//    auto vm = vld1q_u32(m);
//    x = vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(x),vm));
//    std::cout << "x\t" << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << x[3]
//    << "\n";
//    }
//
//    x = vaddq_f32(x,shuf);
//    std::cout << "x\t" << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << x[3]
//    << "\n";
//
//    return x;
//}
// #endif
//
//
////hadamardx8
//// | 1           1          || 1     1                || 1  1 | / |    1 1 ||
/// 1     1             || 1 -1                   | / |       1           1 ||
/// 1    -1                ||       1  1             | / |          1 1
///||    1    -1             ||       1 -1             | / | 1          -1 || 1
/// 1    ||             1  1       | / |    1          -1       || 1     1 || 1
///-1       | / |       1          -1    ||             1    -1    || 1  1 | / |
/// 1          -1 ||                1    -1 ||                   1 -1 |
////
// #if defined(__AVX__)
// inline simd<float,8> hadamard(__m256 x)
//{
//     constexpr auto kMask = floatMask<float>::kSign;
//
//     // x * sqrt(1/8)
//     x = _mm256_mul_ps(x, _mm256_set1_ps(0.3535533905932738f));
//
//     auto shuf = _mm256_shuffle_ps(x,x,_MM_SHUFFLE(2,3,0,1));
//     auto m = _mm256_set_ps(kMask.f,0,kMask.f,0,kMask.f,0,kMask.f,0);
//     x = _mm256_xor_ps(x,m);
//
//     x = _mm256_add_ps(x,shuf);
//
//     shuf = _mm256_shuffle_ps(x,x,_MM_SHUFFLE(1,0,3,2));
//     m = _mm256_set_ps(kMask.f,kMask.f,0,0,kMask.f,kMask.f,0,0);
//     x = _mm256_xor_ps(x,m);
//
//     x = _mm256_add_ps(x,shuf);
//
//     shuf = _mm256_permute2f128_ps(x,x,1);
//     m = _mm256_set_ps(kMask.f,kMask.f,kMask.f,kMask.f,0,0,0,0);
//     x = _mm256_xor_ps(x,m);
//
//     x = _mm256_add_ps(x,shuf);
//
//     return x;
//
//     // is this faster?
//     //auto xl = _mm256_castps256_ps128(x);
//     //auto xh = _mm256_extractf128_ps(x,1);
//
//     //auto a0 = _mm_hadd_ps(xl,xh);
//     //auto a1 = _mm_hsub_ps(xl,xh);
//
//     //auto b0 = _mm_hadd_ps(a0,a1);
//     //auto b1 = _mm_hsub_ps(a0,a1);
//
//     //auto c0 = _mm_hadd_ps(b0,b1);
//     //auto c1 = _mm_hsub_ps(b0,b1);
//
//     //x = _mm256_castps128_ps256(c0);
//     //x = _mm256_insertf128_ps(x,c1,1);
//
//     //return x;
// }
// #endif

// template <typename F, size_t N, size_t H = ilog2(N)>
// Sample<F, N> hadamardBase(Sample<F, N> x)
//{
//     /* compute hadamard transform of sample (without power normalization) */
//     static constexpr auto kL = 1 << H;
//
//     if constexpr (H == 0) {
//         return x;
//     } else {
//         static_assert(N % kL == 0, "N must be divisible by L");
//         constexpr auto kL2 = kL / 2;
//
//         decltype(x) y;
//
// #pragma omp simd
//         for (size_t i = 0; i < N; i += kL) {
//             for (size_t j = 0; j < kL2; ++j) {
//                 auto n     = i + j;
//                 y[n]       = x[n] + x[n + kL2];
//                 y[n + kL2] = x[n] - x[n + kL2];
//             }
//         }
//         //  x0  x0  x2  x2  x4  x4  x6  x6
//         //  x1 -x1  x3 -x3  x5 -x5  x7 -x7
//         //
//         //  x0  x1  x0  x1  x4  x5  x4  x5
//         //  x2  x3 -x2 -x3  x6  x7 -x6 -x7
//         //
//         //  x0  x1  x2  x3  x0  x1  x2  x3
//         //  x4  x5  x6  x7 -x4 -x5 -x6 -x7
//         return hadamardBase<F, N, H - 1>(y);
//     }
// }
//
// template <typename F, size_t N, size_t H = ilog2(N)>
// Sample<F, N> hadamard(Sample<F, N> x)
//{
//     /* hadamard transform with power normalization */
//
//     // first we multiply values by power normalization
//     const F kPowerNorm = std::pow(F(1) / F(2), F(H) / F(2));
//     for (size_t i = 0; i < N; ++i) {
//         x[i] *= kPowerNorm;
//     }
//
//     return hadamardBase<F, N, H>(x);
// }
//
// template <size_t N>
// static inline linalg::fMatrix<N> hadamardInterpolMatrix(float t)
//{
//     /*
//      *  We want to interpolate between identity matrix and hadamard matrix
//      *  ie a matrix H_(t) such as H_(0) = Id, H_(1) = H and H_(t) is
//      orthogonal
//      *  for every t in [0,1]
//      *
//      *  https://dafx2020.mdw.ac.at/proceedings/papers/DAFx2020_paper_53.pdf
//      *  A solution is
//      *
//      *  H_(t) = exp(t⋅log(H))
//      *
//      *  We can diagonalise H = UDU^T, then
//      *  H_(t) = Uexp(t⋅log(D))U^T
//      *
//      *      | 1  1  1  1 |
//      *      | 1 -1  1 -1 |
//      * H=1/2| 1  1 -1 -1 |
//      *      | 1 -1 -1  1 |
//      *
//      * 4 eigen vectors
//      *     v1 = |3  1  1  1| / sqrt(12) eigen value 1
//      *     v2 = |0  1  1 -2| / sqrt(6)  eigen value 1
//      *     v3 = |1 -1 -1 -1| / 2        eigen value -1
//      *     v4 = |0  1 -1  0| / sqrt(2)  eigen value -1
//      *
//      *     |1  0  0  0|
//      * D = |0  1  0  0|  and U = |v1 v2 v3 v4|
//      *     |0  0 -1  0|
//      *     |0  0  0 -1|
//      *
//      *
//      *          |0  0  0   0|                  |1  0    0          0      |
//      * log(D) = |0  0  0   0|  exp(t⋅log(D)) = |0  1    0          0      |
//      *          |0  0  0 -pi|                  |0  0 cos(t⋅pi) -sin(t⋅pi) |
//      *          |0  0 pi   0|                  |0  0 sin(t⋅pi)  cos(t⋅pi) |
//      *
//      * see
//      * https://www.m-hikari.com/ija/ija-password-2008/ija-password1-4-2008/morsyIJA1-4-2008.pdf
//      *
//      * with
//      *                        | 3/sqrt(12)     0      |
//      * A = |v1 v2||v1 v2|^T = | 1/sqrt(12)  1/sqrt(6) | | 3/sqrt(12)
//      1/sqrt(12)
//      * 1/sqrt(12) 1/sqrt(12) | | 1/sqrt(12)  1/sqrt(6) | |     0 1/sqrt(6)
//      * 1/sqrt(6)  -2/sqrt(6) | | 1/sqrt(12) -2/sqrt(6) |
//      *
//      *                        | 3/4  1/4  1/4  1/4 |
//      *                      = | 1/4  1/4  1/4 -1/4 |
//      *                        | 1/4  1/4  1/4 -1/4 |
//      *                        | 1/4 -1/4 -1/4  3/4 |
//      *
//      *
//      *                        |  1/2  0        |
//      * B = |v3 v4||v3 v4|^T = | -1/2  1/sqrt(2)|| 1/2   -1/2       -1/2 -1/2|
//      *                        | -1/2 -1/sqrt(2)|| 0   1/sqrt(2) -1/sqrt(2) 0|
//      *                        | -1/2  0        |
//      *
//      *                        |  1/4 -1/4 -1/4 -1/4 |
//      *                        | -1/4  3/4 -1/4  1/4 |
//      *                        | -1/4 -1/4  3/4  1/4 |
//      *                        | -1/4  1/4  1/4  1/4 |
//      *
//      *                         |  1/2  0        |
//      * C = |v3 v4||-v4 v3|^T = | -1/2  1/sqrt(2)|| 0   -1/sqrt(2)  1/sqrt(2)
//      0|
//      *                         | -1/2 -1/sqrt(2)|| 1/2   -1/2       -1/2
//      -1/2|
//      *                         | -1/2  0        |
//      *
//      *                         |     0      -1/sqrt(8)   1/sqrt(8)   0       
//      |
//      *                         | 1/sqrt(8)     0        -1/sqrt(2)
//      -1/sqrt(8)|
//      *                         |-1/sqrt(8)   1/sqrt(2)     0 1/sqrt(8)|
//      *                         |     0       1/sqrt(8)  -1/sqrt(8)    0 |
//      *
//      *
//      * H_(t) = A + cos(t⋅pi)⋅B + sin(t⋅pi)⋅C
//      *
//      */
//
//     static_assert(N == 4, "Implemented only for N=4");
//     linalg::fMatrix<N> matA = {{{
//         {0.75f, 0.25f, 0.25f, 0.25f},
//         {0.25f, 0.25f, 0.25f, -0.25f},
//         {0.25f, 0.25f, 0.25f, -0.25f},
//         {0.25f, -0.25f, -0.25f, 0.75f},
//     }}};
//     linalg::fMatrix<N> matB = {{{
//         {0.25f, -0.25f, -0.25f, -0.25f},
//         {-0.25f, 0.75f, -0.25f, 0.25f},
//         {-0.25f, -0.25f, 0.75f, 0.25f},
//         {-0.25f, 0.25f, 0.25f, 0.25f},
//     }}};
//
//     constexpr auto kV12     = constants<float>::sqrt1_2;
//     constexpr auto kV18     = 0.5f * kV12;
//     linalg::fMatrix<N> matC = {{{
//         {0.f, kV18, -kV18, 0.f},
//         {-kV18, 0.f, kV12, kV18},
//         {kV18, -kV12, 0.f, -kV18},
//         {0.f, -kV18, kV18, 0.f},
//     }}};
//
//     linalg::fMatrix<N> matH{};
//
//     auto c = cosf(dsp::constants<float>::pi * t);
//     auto s = sinf(dsp::constants<float>::pi * t);
// #pragma omp simd
//     for (size_t j = 0; j < N; ++j) {
//         for (size_t i = 0; i < N; ++i) {
//             matH[j][i] = matA[j][i] + c * matB[j][i] + s * matC[j][i];
//         }
//     }
//
//     return matH;
// }
//
} // namespace dsp
