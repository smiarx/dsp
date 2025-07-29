#pragma once

#include "Context.h"
#include "FastMath.h"
#include <cmath>

namespace dsp
{

template <typename T, int Order> class IIRFilter
{
    using bT = baseType<T>;

  public:
    static constexpr auto kNsos = Order / 2;
    IIRFilter() { static_assert(Order % 2 == 0, "only pair order iirfilter"); }

    class SOS
    {
        /* iir second order section */
      public:
        struct State : public std::array<T, 2> {
            State() : std::array<T, 2>{} {}
        };

        // transfert function type
        using tf = std::array<std::array<bT, 3>, 2>;

        template <class Ctxt, class Stte>
        void process(Ctxt ctxt, Stte &state) const;

        auto tfAnalog(const tf &ba, const T &c, const T &csq);

      private:
        T b_[2];
        T a_[2];
    };

    using State = std::array<typename SOS::State, kNsos>;

    using tf = std::array<typename SOS::tf, kNsos>;

    template <class Ctxt, class State> void process(Ctxt c, State &state) const
    {
        static_assert(!Ctxt::kUseVec, "IIRFilter cannot be vectorized");
        auto x = c.getInput();

        // set gain
        x *= load(b0_);
        c.setOutput(x);

#pragma omp simd
        for (int nsos = 0; nsos < kNsos; ++nsos) {
            sos_[nsos].process(c, state[nsos]);
        }
    }

    void tfAnalog(const tf &ba, const T &freq)
    {
        b0_ = 1;

        /* bilinear transform */
        auto c   = tan(dsp::constants<T>::pi * bT(0.5) * load(freq));
        auto csq = c * c;

        for (int nsos = 0; nsos < kNsos; ++nsos) {
            auto factor = sos_[nsos].tfAnalog(ba[nsos], c, csq);
            b0_ *= load(factor);
        }
    }

  private:
    template <bool highpass> static constexpr void butterworthTF(tf &ba)
    {
        for (int nsos = 0; nsos < kNsos; ++nsos) {
            ba[nsos][0][0] = highpass ? 1 : 0;
            ba[nsos][0][1] = 0;
            ba[nsos][0][2] = highpass ? 0 : 1;
            ba[nsos][1][0] = 1;
            ba[nsos][1][1] =
                -2 * std::cos(dsp::constants<bT>::pi * (2 * nsos + Order + 1) /
                              (2 * Order));
            ba[nsos][1][2] = 1;
        }
    }

  public:
    void butterworthLP(const T &freq)
    {
        tf ba;
        butterworthTF<false>(ba);
        tfAnalog(ba, freq);
    }
    void butterworthHP(const T &freq)
    {
        tf ba;
        butterworthTF<true>(ba);
        tfAnalog(ba, freq);
    }

  private:
    T b0_;
    SOS sos_[kNsos];
};

template <typename T, int Order>
template <class Ctxt, class Stat>
void IIRFilter<T, Order>::SOS::process(Ctxt ctxt, Stat &state) const
{
    /* Transposed Direct Form II */
    auto x = ctxt.getInput();

    decltype(x) s[2];

    auto y = state[0] + x;
    for (int j = 0; j < 2; ++j) s[j] = b_[j] * x - a_[j] * y;

    state[0] = state[1] + s[0];
    state[1] = s[1];

    ctxt.setOutput(y);
}

template <typename T, int Order>
auto IIRFilter<T, Order>::SOS::tfAnalog(const tf &ba, const T &c, const T &csq)
{
    /* filter params */
    /* we define a filter designed with analog coefficients
     * and use bilinear transform to find the corresponding digital coefficients
     * for the desired frequency
     *
     * we use second order section filters (sos) for stability
     */
    T factor;

    assert(ba[1][0] == bT(1));
    auto a0 = ba[1][2];
    auto a1 = ba[1][1];
    auto b0 = ba[0][2];
    auto b1 = ba[0][1];
    auto b2 = ba[0][0];

    auto vc   = load(c);
    auto vcsq = load(csq);

    auto d = bT(1) / (bT(1) + a1 * vc + a0 * vcsq);

    auto b0d    = (b2 + b1 * vc + b0 * vcsq) * d;
    auto invb0d = bT(1) / b0d;
    auto b1d    = bT(2) * (-b2 + b0 * vcsq) * d * invb0d;
    auto b2d    = (b2 - b1 * vc + b0 * vcsq) * d * invb0d;
    auto a1d    = bT(2) * (-bT(1) + a0 * vcsq) * d;
    auto a2d    = (bT(1) - a1 * vc + a0 * vcsq) * d;

    b_[0]  = b1d;
    b_[1]  = b2d;
    a_[0]  = a1d;
    a_[1]  = a2d;
    factor = b0d;

    return factor;
}
} // namespace dsp
