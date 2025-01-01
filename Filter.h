#pragma once

#include "Signal.h"
#include <cassert>
#include <cmath>


namespace dsp
{
template <int N, int NSoS> class Filter;

template <int N> class FilterSoS
{
    template <int Ns, int NSoS> friend class Filter;
    /* second order section */
  public:
    enum class SosPos
    {
        First,
        Rest,
    };

    template <SosPos Pos = SosPos::First> void process(Signal<N> &__restrict x);

  private:
    Signal<N> b[3];
    Signal<N> a[3];
    Signal<N> s[2];
};

template <int N, int NSoS> class Filter
{
  public:
    template <bool NFreq = true>
    static Filter sosanalog(float ba[NSoS][2][3], Signal<NFreq ? N : 1> freq,
                            float samplerate);
    template <bool NFreq = true>
    static Filter butterworthLP(Signal<NFreq ? N : 1>freq, float samplerate);

    void process(Signal<N> &__restrict x)
    {
        sos[0].template process(x);
        for (int j = 1; j < NSoS; ++j) {
            sos[j].template process<FilterSoS<N>::SosPos::Rest>(x);
        }
    }

  private:
    FilterSoS<N> sos[NSoS];
};

template <int N, int NSoS>
template <bool NFreq>
Filter<N, NSoS> Filter<N, NSoS>::sosanalog(float ba[NSoS][2][3],
                                           Signal<NFreq ? N : 1> freq,
                                           float samplerate)
{
    /* filter params */
    /* we define a filter designed with analog coefficients
     * and use bilinear transform to find the corresponding digital coefficients
     * for the desired frequency
     *
     * we use second order section filters (sos) for stability
     */

    constexpr auto NCoef = NFreq ? N : 1;
    using Coef = Signal<NCoef>;

    Filter<N, NSoS> filter;
    for (int j = 0; j < NSoS; ++j) {

        /* bilinear transform */
        Coef c;
        Coef csq;
        for (int i = 0; i < NCoef; ++i) {
            float w1 = 2.f * M_PIf * freq[i];
            c[i]     = 1 / tanf(w1 * 0.5f / samplerate);
            csq[i]   = c[i] * c[i];
        }

        for (int j = 0; j < NSoS; ++j) {
            auto &sos = filter.sos[j];
            for (int i = 0; i < NCoef; ++i) {
                assert(ba[j][1][0] == 1.f);
                float a0  = ba[j][1][2];
                float a1  = ba[j][1][1];
                float b0  = ba[j][0][2];
                float b1  = ba[j][0][1];
                float b2  = ba[j][0][0];
                float d   = 1.f / (a0 + a1 * c[i] + csq[i]);
                float b0d = (b0 + b1 * c[i] + b2 * csq[i]) * d;
                float b1d = 2 * (b0 - b2 * csq[i]) * d;
                float b2d = (b0 - b1 * c[i] + b2 * csq[i]) * d;
                float a1d = 2 * (a0 - csq[i]) * d;
                float a2d = (a0 - a1 * c[i] + csq[i]) * d;

                sos.b[0][i] = b0d;
                sos.b[1][i] = b1d;
                sos.b[2][i] = b2d;
                sos.a[0][i] = 1.f;
                sos.a[1][i] = a1d;
                sos.a[2][i] = a2d;
            }
        }
    }

    /* normalize so that all discret b0 except the first one
     * is equal to 1.0 */
    auto &sos0 = filter.sos[0];
    for (int j = 1; j < NSoS; ++j) {
        auto &sosj = filter.sos[j];
        for (int i = 0; i < NCoef; ++i) {
            auto b0    = sosj.b[0][i];
            auto b0div = 1.f / b0;
            sos0.b[0][i] *= b0;
            sos0.b[1][i] *= b0;
            sos0.b[2][i] *= b0;
            sosj.b[0][i] = 1.f;
            sosj.b[1][i] *= b0div;
            sosj.b[2][i] *= b0div;
        }
    }

    if (!NFreq) {
        /* copy all coefs */
        for (int j = 0; j < NSoS; ++j) {
            auto &sos = filter.sos[j];
            float b0  = sos.b[0][0];
            float b1  = sos.b[1][0];
            float b2  = sos.b[2][0];
            float a0  = sos.a[0][0];
            float a1  = sos.a[1][0];
            float a2  = sos.a[2][0];
            for (int i = 0; i < N; ++i) {
                sos.b[0][i] = b0;
                sos.b[1][i] = b1;
                sos.b[2][i] = b2;
                sos.a[0][i] = a0;
                sos.a[1][i] = a1;
                sos.a[2][i] = a2;
            }
        }
    }
    return filter;
}

template <int N, int NSoS>
template <bool NFreq>
Filter<N, NSoS> Filter<N, NSoS>::butterworthLP(
                                           Signal<NFreq ? N : 1> freq,
                                           float samplerate)
{
    constexpr auto NCoef = NFreq ? N : 1;
    using Coef = Signal<NCoef>;

    Filter<N, NSoS> filter;
    auto& sos = filter.sos[0];

    for(int i = 0; i < NCoef; ++i)
    {
        float c   = tanf(M_PIf * freq[i] / samplerate);
        float csq = c * c;
        float d   = 1.f / (1.f + sqrtf(2.f) * c + csq);

        float a0 = 1.f;
        float a1 = 2.f * (csq - 1.f) * d;
        float a2 = (1.f - sqrtf(2.f) * c + csq) * d;
        float b0 = csq * d;
        float b1 = 2.f * b0;
        float b2 = b0;

        sos.a[0][i] = a0;
        sos.a[1][i] = a1;
        sos.a[2][i] = a2;
        sos.b[0][i] = b0;
        sos.b[1][i] = b1;
        sos.b[2][i] = b2;
    }

    if (!NFreq) {
    /* copy all coefs */
        float b0  = sos.b[0][0];
        float b1  = sos.b[1][0];
        float b2  = sos.b[2][0];
        float a0  = sos.a[0][0];
        float a1  = sos.a[1][0];
        float a2  = sos.a[2][0];
        for (int i = 0; i < N; ++i) {
            sos.b[0][i] = b0;
            sos.b[1][i] = b1;
            sos.b[2][i] = b2;
            sos.a[0][i] = a0;
            sos.a[1][i] = a1;
            sos.a[2][i] = a2;
        }
    }
    return filter;
}

template <int N>
template <typename FilterSoS<N>::SosPos Pos>
void FilterSoS<N>::process(Signal<N> &__restrict x)
{
/* use direct form II */
#pragma omp simd
    for (int i = 0; i < N; ++i) {
        float a1 = a[1][i];
        float a2 = a[2][i];
        float b0 = b[0][i];
        float b1 = b[1][i];
        float b2 = b[2][i];

        float s1 = s[0][i];
        float s2 = s[1][i];
        float s0 = x[i] - a1 * s1 - a2 * s2;

        if (Pos == SosPos::First) x[i] = b0 * s0 + b1 * s1 + b2 * s2;
        else
            x[i] = s0 + b1 * s1 + b2 * s2;

        s[1][i] = s[0][i];
        s[0][i] = s0;
    }
}

//        /*
// #ifdef FILTER_VECSIZE
//
// typedef struct {
//    float b[3][FILTER_VECSIZE]
//        __attribute__((aligned(sizeof(float) * FILTER_VECSIZE)));
//    float a[3][FILTER_VECSIZE]
//        __attribute__((aligned(sizeof(float) * FILTER_VECSIZE)));
//    float s[2][FILTER_VECSIZE]
//        __attribute__((aligned(sizeof(float) * FILTER_VECSIZE)));
//} filter(_t);
//
// void filter(_butterworth_lp)(filter(_t) * filter, float samplerate, float
// fr); void filter(_sos_analog)(filter(_t) filter[], const float ba[][2][3],
//                         const float freqs[FILTER_VECSIZE], float samplerate,
//                         int nsos);
//
// #ifdef _FILTER_C
// void filter(_butterworth_lp)(filter(_t) * filter, float samplerate, float fr)
//{
//    float c   = tanf(M_PI * fr / samplerate);
//    float csq = c * c;
//    float d   = 1.f / (1.f + sqrtf(2.f) * c + csq);
//
//    float a0 = 1.f;
//    float a1 = 2.f * (csq - 1.f) * d;
//    float a2 = (1.f - sqrtf(2.f) * c + csq) * d;
//    float b0 = csq * d;
//    float b1 = 2.f * b0;
//    float b2 = b0;
//    for (int i = 0; i < FILTER_VECSIZE; ++i)
//        filter->a[0][i] = a0, filter->a[1][i] = a1, filter->a[2][i] = a2,
//        filter->b[0][i] = b0, filter->b[1][i] = b1, filter->b[2][i] = b2;
//}
//
// void filter(_sos_analog)(filter(_t) filter[], const float ba[][2][3],
//                         const float freqs[FILTER_VECSIZE], float samplerate,
//                         int nsos)
//{
//    /* filter params */
//    /* we define a filter designed with analog coefficients
//     * and use bilinear transform to find the corresponding digital
//     coefficients
//     * for the desired frequency
//     *
//     * we use second order section filters (sos) for stability
//     */
//
//    float c[FILTER_VECSIZE];
//    float csq[FILTER_VECSIZE];
//    /* bilinear transform */
//    for (int i = 0; i < FILTER_VECSIZE; ++i) {
//        float w1 = 2.f * M_PI * freqs[i];
//        c[i]     = 1 / tanf(w1 * 0.5f / samplerate);
//        csq[i]   = c[i] * c[i];
//    }
//    for (int j = 0; j < nsos; ++j) {
//        for (int i = 0; i < FILTER_VECSIZE; ++i) {
//            assert(ba[j][1][0] == 1.f);
//            float a0  = ba[j][1][2];
//            float a1  = ba[j][1][1];
//            float b0  = ba[j][0][2];
//            float b1  = ba[j][0][1];
//            float b2  = ba[j][0][0];
//            float d   = 1.f / (a0 + a1 * c[i] + csq[i]);
//            float b0d = (b0 + b1 * c[i] + b2 * csq[i]) * d;
//            float b1d = 2 * (b0 - b2 * csq[i]) * d;
//            float b2d = (b0 - b1 * c[i] + b2 * csq[i]) * d;
//            float a1d = 2 * (a0 - csq[i]) * d;
//            float a2d = (a0 - a1 * c[i] + csq[i]) * d;
//
//            filter[j].b[0][i] = b0d;
//            filter[j].b[1][i] = b1d;
//            filter[j].b[2][i] = b2d;
//            filter[j].a[0][i] = 1.f;
//            filter[j].a[1][i] = a1d;
//            filter[j].a[2][i] = a2d;
//        }
//    }
//
//    /* normalize so that all discret b0 except the first one
//     * is equal to 1.0 */
//    for (int j = 1; j < nsos; ++j) {
//        for (int i = 0; i < FILTER_VECSIZE; ++i) {
//            float b0 = filter[j].b[0][i];
//            filter[0].b[0][i] *= b0;
//            filter[0].b[1][i] *= b0;
//            filter[0].b[2][i] *= b0;
//            float b0div       = 1.f / b0;
//            filter[j].b[0][i] = 1.f;
//            filter[j].b[1][i] *= b0div;
//            filter[j].b[2][i] *= b0div;
//        }
//    }
//}
// #endif
//
// static inline void filter(_process)(filter(_t) * filter,
//                                    float y[restrict FILTER_VECSIZE], int
//                                    nsos)
//{
//    y = (float *)__builtin_assume_aligned(y, sizeof(float) * FILTER_VECSIZE);
//
//    for (int j = 0; j < nsos; ++j) {
///* use direct form II */
// #pragma omp simd
//         for (int i = 0; i < FILTER_VECSIZE; ++i) {
//             float a1 = filter[j].a[1][i];
//             float a2 = filter[j].a[2][i];
//             float b0 = filter[j].b[0][i];
//             float b1 = filter[j].b[1][i];
//             float b2 = filter[j].b[2][i];
//
//             float s1 = filter[j].s[0][i];
//             float s2 = filter[j].s[1][i];
//             float s0 = y[i] - a1 * s1 - a2 * s2;
//
//             if (j == 0) y[i] = b0 * s0 + b1 * s1 + b2 * s2;
//             else
//                 y[i] = s0 + b1 * s1 + b2 * s2;
//
//             filter[j].s[1][i] = filter[j].s[0][i];
//             filter[j].s[0][i] = s0;
//         }
//     }
// }
// #endif
//*/
}
