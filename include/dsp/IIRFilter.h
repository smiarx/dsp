#pragma once

#include "Context.h"
#include "FastMath.h"
#include "Utils.h"
#include <cmath>

#include "gcem/gcem.h"

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

// second order section
template <typename T> struct SOS {
    T b2, b1, b0, a2, a1, a0;
};

template <typename T, size_t K> class NormalizedBiquadSeries
{
    /* biquad iir filter in series
     * all biquads have b0 == 1
     * */
    using bT                  = baseType<T>;
    static constexpr auto kN  = kTypeWidth<T>;
    static constexpr auto kNK = kN * K;

    using simd1 = simd<bT, kNK>;

  public:
    struct State : public std::array<multi<T, K>, 2> {
        State() : std::array<multi<T, K>, 2>{} {}
    };

    // second order section type
    struct SOS {
        multi<T, K> b2{}, b1{}, b0{}, a2{}, a1{}, a0{};
    };

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

    auto fromAnalog(const SOS sos, const T &c, const T &csq)
    {
        /* filter params */
        /* we define a filter designed with analog coefficients
         * and use bilinear transform to find the corresponding digital
         * coefficients for the desired frequency
         *
         * we use second order section filters (sos) for stability
         */
        auto a1 = load(sos.a1);
        auto a0 = load(sos.a0);
        auto b2 = load(sos.b2);
        auto b1 = load(sos.b1);
        auto b0 = load(sos.b0);
        assert(all(load(sos.a2).operator==(bT(1))));

        auto d = bT(1) / (bT(1) + a1 * c + a0 * csq);

        auto b0d    = (b2 + b1 * c + b0 * csq) * d;
        auto invb0d = bT(1) / b0d;

        auto b1d = bT(2) * (-b2 + b0 * csq) * d * invb0d;
        auto b2d = (b2 - b1 * c + b0 * csq) * d * invb0d;
        auto a1d = bT(2) * (-bT(1) + a0 * csq) * d;
        auto a2d = (bT(1) - a1 * c + a0 * csq) * d;

        b_[0] = b1d;
        b_[1] = b2d;
        a_[0] = a1d;
        a_[1] = a2d;

        return product<kN>(b0d);
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

    // biquad type
    using Biquad = NormalizedBiquadSeries<T, kNsos>;

  public:
    static constexpr auto kTotalNsos = Order / 2;
    IIRFilter() { static_assert(Order % 2 == 0, "only pair order iirfilter"); }

    // second order section type
    // similar to scipy.signal sos type
    struct SOS {
        constexpr SOS(std::array<dsp::SOS<bT>, Order / 2> sos) : SOS(sos.data())
        {
        }
        constexpr SOS(dsp::SOS<bT> *sos) : sosNext(sos + kThisOrder / 2)
        {
            for (size_t k = 0; k < kThisOrder / 2; ++k) {
                sos0.b2[k] = sos[k].b2;
                sos0.b1[k] = sos[k].b1;
                sos0.b0[k] = sos[k].b0;
                sos0.a2[k] = sos[k].a2;
                sos0.a1[k] = sos[k].a1;
                sos0.a0[k] = sos[k].a0;
            }
            for (size_t k = kThisOrder / 2; k < kNsos; ++k) {
                sos0.a2[k] = 1;
                sos0.b2[k] = 1;
            }
        }
        typename Biquad::SOS sos0{};
        typename NextIIR::SOS sosNext;
    };

    struct State {
        typename NormalizedBiquadSeries<T, kNsos>::State s0;
        typename NextIIR::State sNext;
    };

    template <class Ctxt> void processBlock(Ctxt c, State &state) const
    {
        static_assert(!Ctxt::kUseVec, "IIRFilter cannot be vectorized");

        if constexpr (!Normalized) {
            CTXTRUNVEC(c)
            {
                auto x = c.getInput();
                x *= gain_;
                c.setOutput(x);
            };
        }

        CTXTRUN(c) { biquad_.process(c, state.s0); };

        if constexpr (Order > kNsos * 2) {
            NextIIR::processBlock(c, state.sNext);
        }
    }

    template <class Ctxt> void process(Ctxt c, State &state) const
    {
        // compute one element with block of size = 1
        c.setBlockSize(1);
        processBlock(c, state);
    }

    // convert from analog scipy sos filter type
    void fromAnalog(const SOS &sos, const T &freq)
    {
        static_assert(!Normalized);

        /* bilinear transform */
        auto c   = tan(dsp::constants<T>::pi * bT(0.5) * load(freq));
        auto csq = c * c;

        auto factor = biquad_.fromAnalog(sos.sos0, c, csq);

        if constexpr (Order > kThisOrder) {
            factor *= NextIIR::fromAnalog(sos.sosNext, c, csq);
        }

        if constexpr (!Normalized) {
            gain_ = factor;
        } else {
            assert(all(factor == T(1)));
        }
    }

  protected:
    auto fromAnalog(const SOS &sos, const T &c, const T &csq)
    {
        static_assert(Normalized);
        auto factor = biquad_.fromAnalog(sos.sos0, c, csq);
        auto b0     = load(factor);
        if constexpr (Order > kNsos * 2) {
            b0 *= NextIIR::fromAnalog(sos.sosNext, c, csq);
        }
        return b0;
    }

  private:
    template <bool highpass> static constexpr SOS butterworthTF()
    {
        std::array<dsp::SOS<bT>, Order / 2> tf{};
        for (size_t nsos = 0; nsos < Order / 2; ++nsos) {
            tf[nsos].b2 = highpass ? 1 : 0;
            tf[nsos].b1 = 0;
            tf[nsos].b0 = highpass ? 0 : 1;
            tf[nsos].a2 = 1;
            tf[nsos].a1 =
                -2 * gcem::cos(dsp::constants<bT>::pi *
                               bT(2 * nsos + Order + 1) / bT(2 * Order));
            tf[nsos].a0 = 1;
        }
        return tf;
    }

  public:
    void butterworthLP(const T &freq)
    {
        constexpr SOS kSos = butterworthTF<false>();
        fromAnalog(kSos, freq);
    }
    void butterworthHP(const T &freq)
    {
        constexpr SOS kSos = butterworthTF<true>();
        fromAnalog(kSos, freq);
    }

  private:
    struct Empty {
    };
    std::conditional_t<Normalized, Empty, T> gain_;
    Biquad biquad_;
};

template <typename T> class IIRFilter<T, 0, true>
{
  protected:
    struct SOS {
        SOS() = default;
        constexpr SOS(dsp::SOS<baseType<T>> *) {}
    };
    struct State {
    };
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
