#pragma once

#include "Context.h"
#include "FastMath.h"
#include "Signal.h"
#include <cmath>

namespace dsp
{

template <int N, int Order> class IIRFilter
{
  public:
    static constexpr auto NSOS = Order / 2;
    IIRFilter() { static_assert(Order % 2 == 0, "only pair order iirfilter"); }

    class SOS
    {
        /* iir second order section */
      public:
        template <int N_ = N> using Mem = fData<N_>[2];

        template <class Ctxt, class Mem>
        void process(Ctxt ctxt, Mem &mem) const;

        fData<N> tfAnalog(const float ba[2][3], fData<N> c, fData<N> csq);

      private:
        fData<N> b_[2];
        fData<N> a_[2];
    };

    template <int N_ = N>
    using Mem = std::array<typename SOS::template Mem<N_>, NSOS>;

    template <class Ctxt, class Mem> void process(Ctxt c, Mem &mem) const
    {
        auto &x = c.getSignal();

#pragma omp simd
        inFor(x, k, i) { x[k][i] *= b0_[i % N]; }
        for (int nsos = 0; nsos < NSOS; ++nsos) {
            sos_[nsos].process(c, mem[nsos]);
        }
    }

    void tfAnalog(const float ba[NSOS][2][3], fData<N> freq)
    {
        arrayFor(b0_, i) { b0_[i] = 1.f; }

        /* bilinear transform */
        fData<N> c;
        fData<N> csq;
        arrayFor(c, i)
        {
            c[i]   = tanf(dsp::constants<float>::pi * 0.5f * freq[i]);
            csq[i] = c[i] * c[i];
        }

        for (int nsos = 0; nsos < NSOS; ++nsos) {
            auto factor = sos_[nsos].tfAnalog(ba[nsos], c, csq);
            arrayFor(b0_, i) { b0_[i] *= factor[i]; }
        }
    }

    template <bool highpass>
    static constexpr void butterworthTF(float ba[NSOS][2][3])
    {
        for (int nsos = 0; nsos < NSOS; ++nsos) {
            ba[nsos][0][0] = highpass ? 1.f : 0.f;
            ba[nsos][0][1] = 0.f;
            ba[nsos][0][2] = highpass ? 0.f : 1.f;
            ba[nsos][1][0] = 1.f;
            ba[nsos][1][1] = -2.f * cosf(dsp::constants<float>::pi *
                                         (2 * nsos + Order + 1) / (2 * Order));
            ba[nsos][1][2] = 1.f;
        }
    }

    void butterworthLP(fData<N> freq)
    {
        float ba[NSOS][2][3];
        butterworthTF<false>(ba);
        tfAnalog(ba, freq);
    }
    void butterworthHP(fData<N> freq)
    {
        float ba[NSOS][2][3];
        butterworthTF<true>(ba);
        tfAnalog(ba, freq);
    }

  private:
    fData<N> b0_;
    SOS sos_[NSOS];
};

template <int N, int Order>
template <class Ctxt, class Memory>
void IIRFilter<N, Order>::SOS::process(Ctxt ctxt, Memory &mem) const
{
    /* Transposed Direct Form II */
    auto &in = ctxt.getSignal();

    arrayFor(in, k)
    {
        /* we double input and ouput to use simd */
        typename Ctxt::BaseType x;
        typename Ctxt::BaseType y;
        typename Ctxt::BaseType s[2];

        x = in[k];

        arrayFor(y, i) { y[i] = mem[0][i] + x[i]; }

        for (int j = 0; j < 2; ++j) {
            arrayFor(y, i)
            {
                s[j][i] = b_[j][i % N] * x[i];
                s[j][i] -= a_[j][i % N] * y[i];
            }
        }
        arrayFor(y, i) { mem[0][i] = mem[1][i] + s[0][i]; }
        arrayFor(y, i) { mem[1][i] = s[1][i]; }

        in[k] = y;
    }
}

template <int N, int Order>
fData<N> IIRFilter<N, Order>::SOS::tfAnalog(const float ba[2][3], fData<N> c,
                                            fData<N> csq)
{
    /* filter params */
    /* we define a filter designed with analog coefficients
     * and use bilinear transform to find the corresponding digital coefficients
     * for the desired frequency
     *
     * we use second order section filters (sos) for stability
     */
    fData<N> factor;

    assert(ba[1][0] == 1.f);
    float a0 = ba[1][2];
    float a1 = ba[1][1];
    float b0 = ba[0][2];
    float b1 = ba[0][1];
    float b2 = ba[0][0];

    for (int i = 0; i < N; ++i) {
        float d = 1.f / (1.f + a1 * c[i] + a0 * csq[i]);

        float b0d    = (b2 + b1 * c[i] + b0 * csq[i]) * d;
        float invb0d = 1.f / b0d;
        float b1d    = 2 * (-b2 + b0 * csq[i]) * d * invb0d;
        float b2d    = (b2 - b1 * c[i] + b0 * csq[i]) * d * invb0d;
        float a1d    = 2 * (-1.f + a0 * csq[i]) * d;
        float a2d    = (1.f - a1 * c[i] + a0 * csq[i]) * d;

        b_[0][i] = b1d;
        b_[1][i] = b2d;
        a_[0][i] = a1d;
        a_[1][i] = a2d;

        factor[i] = b0d;
    }
    return factor;
}
} // namespace dsp
