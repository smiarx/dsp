#pragma once

#include "Context.h"
#include "FastMath.h"
#include "Utils.h"
#include <cmath>

namespace dsp
{

template <typename T, size_t K> class NormalizedBiquadSeries
{
    /* biquad iir filter in series
     * all biquads have b0 == 1
     * */
    using bT                  = baseType<T>;
    static constexpr auto kN  = kTypeWidth<T>;
    static constexpr auto kNK = kN * K;

  public:
    struct State : public std::array<multi<T, K>, 2> {
        State() : std::array<multi<T, K>, 2>{} {}
    };

    // transfert function type
    using tf = std::array<std::array<bT, 3>, 2>;

    template <class Ctxt> void process(Ctxt ctxt, State &state) const
    {
        /* Transposed Direct Form II */
        auto in = ctxt.getInput();

        auto y = simd<bT, kNK>::convert(in) + prefix<kN>(load(state[0]));
        auto x = push(y, in);

        if constexpr (K < DSP_MAX_VEC_SIZE / sizeof(T)) {
            using simd2 = simd<bT, kNK * 2>;
            auto x2     = simd2::convert(x);
            auto y2     = simd2::convert(y);
            auto b      = simd2::loadu((bT *)&b_[0]);
            auto a      = simd2::loadu((bT *)&a_[0]);
            auto st     = simd2::loadu((bT *)&state[0]);
            st          = shift<-static_cast<int>(kNK)>(st);
            st += b * x2 - a * y2;
            st.storeu((bT *)&state[0]);
        } else {
            auto s0 = b_[0] * x - a_[0] * y;
            auto s1 = b_[1] * x - a_[1] * y;

            state[0] = state[1] + s0;
            state[1] = s1;
        }

        ctxt.setOutput(getlane<K - 1, kN>(y));
    }

    template <size_t Nsos>
    auto tfAnalog(const tf *bacoeffs, const T &c, const T &csq)
    {
        /* filter params */
        /* we define a filter designed with analog coefficients
         * and use bilinear transform to find the corresponding digital
         * coefficients for the desired frequency
         *
         * we use second order section filters (sos) for stability
         */
        auto factor = T(1);

        for (size_t k = 0; k < Nsos; ++k) {
            auto &ba = *bacoeffs;
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

            b_[0][k] = b1d;
            b_[1][k] = b2d;
            a_[0][k] = a1d;
            a_[1][k] = a2d;
            factor *= b0d;

            ++bacoeffs;
        }

        return factor;
    }

  private:
    multi<T, K> b_[2]{};
    multi<T, K> a_[2]{};
};

template <typename T> static constexpr auto getBiquadSizeFromOrder(size_t order)
{
    order /= 2;
    constexpr auto kMaxK = DSP_MAX_VEC_SIZE / sizeof(T);
    if (order >= kMaxK) return kMaxK;
    return nextPow2(order);
}

template <typename T, size_t Order, bool Normalized = false>
class IIRFilter
    : public IIRFilter<
          T, Order - std::min(Order, 2 * getBiquadSizeFromOrder<T>(Order)),
          true>
{
    static_assert(Order % 2 == 0, "Order has to be a multiple of 2");

    using bT = baseType<T>;

    static constexpr auto kNsos      = getBiquadSizeFromOrder<T>(Order);
    static constexpr auto kThisOrder = std::min(Order, 2 * kNsos);
    using NextIIR                    = IIRFilter<T, Order - kThisOrder, true>;

  public:
    static constexpr auto kTotalNsos = Order / 2;
    IIRFilter() { static_assert(Order % 2 == 0, "only pair order iirfilter"); }

    // transfert function type
    using tf1 = std::array<std::array<bT, 3>, 2>;
    using tf  = std::array<tf1, kTotalNsos>;

    struct State {
        typename NormalizedBiquadSeries<T, kNsos>::State s0;
        typename NextIIR::State sNext;
    };

    template <class Ctxt> void process(Ctxt c, State &state) const
    {
        static_assert(!Ctxt::kUseVec, "IIRFilter cannot be vectorized");

        if constexpr (!Normalized) {
            auto x = c.getInput();

            // set gain
            x *= load(b0_);
            c.setOutput(x);
        }

        biquad_.process(c, state.s0);
        if constexpr (Order > kNsos * 2) {
            NextIIR::process(c, state.sNext);
        }
    }

    void tfAnalog(const tf &ba, const T &freq)
    {
        static_assert(!Normalized);

        /* bilinear transform */
        auto c   = tan(dsp::constants<T>::pi * bT(0.5) * load(freq));
        auto csq = c * c;

        auto factor =
            biquad_.template tfAnalog<kThisOrder / 2>(ba.data(), c, csq);

        if constexpr (Order > kThisOrder) {
            factor *= NextIIR::tfAnalog(&ba[kNsos], c, csq);
        }

        if constexpr (!Normalized) {
            b0_ = factor;
        } else {
            assert(all(factor == T(1)));
        }
    }

  protected:
    auto tfAnalog(const tf1 *ba, const T &c, const T &csq)
    {
        static_assert(Normalized);
        auto factor = biquad_.template tfAnalog<kThisOrder / 2>(ba, c, csq);
        auto b0     = load(factor);
        if constexpr (Order > kNsos * 2) {
            b0 *= NextIIR::tfAnalog(&ba[kNsos], c, csq);
        }
        return b0;
    }

  private:
    template <bool highpass> static constexpr void butterworthTF(tf &ba)
    {
        for (size_t nsos = 0; nsos < Order / 2; ++nsos) {
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
    struct Empty {
    };
    std::conditional_t<Normalized, Empty, T> b0_;
    NormalizedBiquadSeries<T, kNsos> biquad_;
};

template <typename T> class IIRFilter<T, 0, true>
{
  protected:
    struct State {
    };
};

} // namespace dsp
