#pragma once

#include "Delay.h"
#include "Signal.h"
#include <cassert>
#include <cmath>

namespace dsp
{
template <int N, int NSoS> class Filter;

template <int N> class FilterSoS
{
    /* second order section */
    template <int Ns, int NSoS> friend class Filter;

  public:
    using DL = CopyDelayLine<N, 2>;
    enum class SosPos {
        First,
        Rest,
    };

    template <class Ctxt, SosPos Pos = SosPos::First>
    void process(Ctxt c, DL &delayline) const
    {
        auto &x = c.getIn();

        /* read delayline */
        auto s1 = TapFix<1>().read(c, delayline);
        auto s2 = TapFix<2>().read(c, delayline);
        decltype(s1) s0;

        /* use direct form II */
        for (int i = 0; i < N; ++i) {
            float a1 = a[1][i];
            float a2 = a[2][i];
            float b0 = b[0][i];
            float b1 = b[1][i];
            float b2 = b[2][i];

            s0[i] = x[i] - a1 * s1[i] - a2 * s2[i];

            if (Pos == SosPos::First)
                x[i] = b0 * s0[i] + b1 * s1[i] + b2 * s2[i];
            else
                x[i] = s0[i] + b1 * s1[i] + b2 * s2[i];
        }
        delayline.write(c, s0);
    }

  private:
    Signal<N> b[3];
    Signal<N> a[3];
};

template <int N, int NSoS = 1> class Filter
{
  public:
    using DL = std::array<typename FilterSoS<N>::DL, NSoS>;

    template <int NFreq>
    void sosanalog(const float ba[NSoS][2][3], Signal<NFreq> freq);
    template <int NFreq> void butterworthLP(Signal<NFreq> freq);

    template <class Ctxt> void process(Ctxt c, DL &delayline) const
    {
        sos[0].process(c, delayline[0]);
        for (int j = 1; j < NSoS; ++j) {
            sos[j].template process<Ctxt, FilterSoS<N>::SosPos::Rest>(
                c, delayline[j]);
        }
    }

  private:
    FilterSoS<N> sos[NSoS];
};

template <int N, int NSoS>
template <int NFreq>
void Filter<N, NSoS>::sosanalog(const float ba[NSoS][2][3], Signal<NFreq> freq)
{
    /* filter params */
    /* we define a filter designed with analog coefficients
     * and use bilinear transform to find the corresponding digital coefficients
     * for the desired frequency
     *
     * we use second order section filters (sos) for stability
     */

    static_assert(NFreq == N || NFreq == 1,
                  "The number of frequency is either the vector size or 1");

    for (int j = 0; j < NSoS; ++j) {

        /* bilinear transform */
        Signal<NFreq> c;
        Signal<NFreq> csq;
        for (int i = 0; i < NFreq; ++i) {
            c[i]   = 1 / tanf(M_PIf * 0.5f * freq[i]);
            csq[i] = c[i] * c[i];
        }

        for (int j = 0; j < NSoS; ++j) {
            for (int i = 0; i < NFreq; ++i) {
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

                sos[j].b[0][i] = b0d;
                sos[j].b[1][i] = b1d;
                sos[j].b[2][i] = b2d;
                sos[j].a[0][i] = 1.f;
                sos[j].a[1][i] = a1d;
                sos[j].a[2][i] = a2d;
            }
        }
    }

    /* normalize so that all discret b0 except the first one
     * is equal to 1.0 */
    auto &sos0 = sos[0];
    for (int j = 1; j < NSoS; ++j) {
        auto &sosj = sos[j];
        for (int i = 0; i < NFreq; ++i) {
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

    /* copy all coefs if NFreq < N */
    for (int j = 0; j < NSoS; ++j) {
        float b0 = sos[j].b[0][0];
        float b1 = sos[j].b[1][0];
        float b2 = sos[j].b[2][0];
        float a0 = sos[j].a[0][0];
        float a1 = sos[j].a[1][0];
        float a2 = sos[j].a[2][0];
        for (int i = NFreq; i < N; ++i) {
            sos[j].b[0][i] = b0;
            sos[j].b[1][i] = b1;
            sos[j].b[2][i] = b2;
            sos[j].a[0][i] = a0;
            sos[j].a[1][i] = a1;
            sos[j].a[2][i] = a2;
        }
    }
}

template <int N, int NSoS>
template <int NFreq>
void Filter<N, NSoS>::butterworthLP(Signal<NFreq> freq)
{
    static_assert(NSoS == 1);
    static_assert(NFreq == N || NFreq == 1,
                  "The number of frequency is either the vector size or 1");

    auto &sos0 = sos[0];

    for (int i = 0; i < NFreq; ++i) {
        float c   = tanf(M_PIf * 0.5f * freq[i]);
        float csq = c * c;
        float d   = 1.f / (1.f + sqrtf(2.f) * c + csq);

        float a0 = 1.f;
        float a1 = 2.f * (csq - 1.f) * d;
        float a2 = (1.f - sqrtf(2.f) * c + csq) * d;
        float b0 = csq * d;
        float b1 = 2.f * b0;
        float b2 = b0;

        sos0.a[0][i] = a0;
        sos0.a[1][i] = a1;
        sos0.a[2][i] = a2;
        sos0.b[0][i] = b0;
        sos0.b[1][i] = b1;
        sos0.b[2][i] = b2;
    }

    /* copy all coefs if NFreq < N*/
    float b0 = sos0.b[0][0];
    float b1 = sos0.b[1][0];
    float b2 = sos0.b[2][0];
    float a0 = sos0.a[0][0];
    float a1 = sos0.a[1][0];
    float a2 = sos0.a[2][0];
    for (int i = NFreq; i < N; ++i) {
        sos0.b[0][i] = b0;
        sos0.b[1][i] = b1;
        sos0.b[2][i] = b2;
        sos0.a[0][i] = a0;
        sos0.a[1][i] = a1;
        sos0.a[2][i] = a2;
    }
}
} // namespace dsp
