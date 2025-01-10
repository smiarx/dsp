#pragma once

#include "Delay.h"
#include "Signal.h"
#include <cassert>
#include <cmath>

namespace dsp
{
enum class SoSPos {
    First,
    Rest,
};

template <int N, int NSoS = 1> class IIRFilter
{
    class SoS
    {
        /* second order section */
        friend class IIRFilter;

      public:
        constexpr SoS()                = default;
        template <int N_ = N> using DL = CopyDelayLine<N_, 2>;

        template <SoSPos Pos = SoSPos::First, class Ctxt, class DL>
        void process(Ctxt c, DL &delayline) const
        {
            auto &x = c.getIn();

            /* read delayline */
            auto s1 = TapFix<1>().read(c, delayline);
            auto s2 = TapFix<2>().read(c, delayline);
            decltype(s1) s0;

            /* use direct form II */
            for (int k = 0; k < x.size(); ++k) {
                for (int i = 0; i < x[0].size(); ++i) {
                    float a1 = a[1][i % N];
                    float a2 = a[2][i % N];
                    s0[k][i] = x[k][i] - a1 * s1[k][i] - a2 * s2[k][i];
                }
            }
            for (int k = 0; k < x.size(); ++k) {
                for (int i = 0; i < x[0].size(); ++i) {
                    float b0 = b[0][i % N];
                    float b1 = b[1][i % N];
                    float b2 = b[2][i % N];

                    if (Pos == SoSPos::First)
                        x[k][i] = b0 * s0[k][i] + b1 * s1[k][i] + b2 * s2[k][i];
                    else
                        x[k][i] = s0[k][i] + b1 * s1[k][i] + b2 * s2[k][i];
                }
            }
            delayline.write(c, s0);
        }

      private:
        Signal<N> b[3] = {0};
        Signal<N> a[3] = {0};
    };

  public:
    template <int N_>
    using DL = std::array<typename SoS::template DL<N_>, NSoS>;

    constexpr IIRFilter()            = default;
    constexpr IIRFilter(IIRFilter &) = default;
    static IIRFilter newButterworthLP(Signal<N> freq)
    {
        IIRFilter f;
        f.butterworthLP(freq);
        return f;
    }

    void sosanalog(const float ba[NSoS][2][3], Signal<N> freq);
    constexpr void butterworthLP(Signal<N> freq);

    template <class Ctxt, class DL> void process(Ctxt c, DL &delayline) const
    {
        sos_[0].process(c, delayline[0]);
        for (int j = 1; j < NSoS; ++j) {
            sos_[j].template process<SoSPos::Rest>(c, delayline[j]);
        }
    }

    // private:
    SoS sos_[NSoS];
};

template <int N, int NSoS>
void IIRFilter<N, NSoS>::sosanalog(const float ba[NSoS][2][3], Signal<N> freq)
{
    /* filter params */
    /* we define a filter designed with analog coefficients
     * and use bilinear transform to find the corresponding digital coefficients
     * for the desired frequency
     *
     * we use second order section filters (sos) for stability
     */

    for (int j = 0; j < NSoS; ++j) {

        /* bilinear transform */
        Signal<N> c;
        Signal<N> csq;
        for (int i = 0; i < N; ++i) {
            c[i]   = 1 / tanf(M_PIf * 0.5f * freq[i]);
            csq[i] = c[i] * c[i];
        }

        for (int j = 0; j < NSoS; ++j) {
            for (int i = 0; i < N; ++i) {
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

                sos_[j].b[0][i] = b0d;
                sos_[j].b[1][i] = b1d;
                sos_[j].b[2][i] = b2d;
                sos_[j].a[0][i] = 1.f;
                sos_[j].a[1][i] = a1d;
                sos_[j].a[2][i] = a2d;
            }
        }
    }

    /* normalize so that all discret b0 except the first one
     * is equal to 1.0 */
    auto &sos0 = sos_[0];
    for (int j = 1; j < NSoS; ++j) {
        auto &sosj = sos_[j];
        for (int i = 0; i < N; ++i) {
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
}

template <int N, int NSoS>
constexpr void IIRFilter<N, NSoS>::butterworthLP(Signal<N> freq)
{
    static_assert(NSoS == 1);

    auto &sos0 = sos_[0];

    for (int i = 0; i < N; ++i) {
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
}
} // namespace dsp
