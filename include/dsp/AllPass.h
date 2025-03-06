#pragma once

#include "Delay.h"
#include "Signal.h"
#include <cmath>

namespace dsp
{
template <class... T> class List
{
};

template <int N, class Tap = TapTail> class AllPass
{
    /* Implement a Allpass filter with coeff a
     * template:
     *      N: vector size
     *      Tap: a Tap class for the recursive delayline
     *
     *  Allpass process needs a DelayLineType delayline to process
     *  internal state
     */
  public:
    constexpr AllPass() = default;
    template <class... Args>
    constexpr AllPass(const fData<N> &a, Args... args) : a_(a), tap_(args...)
    {
    }

    template <class Ctxt, class DL> void process(Ctxt c, DL &delayline) const
    {
        static_assert(Ctxt::VecSize <= DL::Length);

        auto &x = c.getSignal();
        auto sN = tap_.read(c, delayline);
        decltype(sN) s0;

        for (size_t k = 0; k < x.size(); ++k) {
            for (size_t i = 0; i < x[0].size(); ++i) {
                auto a = a_[i % N];

                s0[k][i] = x[k][i] - a * sN[k][i];
                x[k][i]  = a * s0[k][i] + sN[k][i];
            }
        }
        delayline.write(c, s0);
    }

    void setCoeff(fData<N> a) { a_ = a; }
    const auto &getCoeff() const { return a_; }
    void setDelay(fData<N> d) { tap_.setDelay(d); }
    void setDelay(iData<N> d) { tap_.setDelay(d); }

    const auto &getTap() const { return tap_; }

  private:
    fData<N> a_;
    Tap tap_;
};

template <int N, class TapOut = TapTail, class TapIn = TapTail>
class TapAllPass : public TapOut
{
  public:
    static constexpr auto minfdelay = 0.5f;

    constexpr TapAllPass() = default;
    template <class... Args>
    constexpr TapAllPass(const fData<N> &a, Args... args) :
        TapOut(args...), allpass_(a)
    {
    }

    void setDelay(fData<N> d)
    {
        fData<N> a;
        iData<N> id;
        for (int i = 0; i < N; ++i) {
            id[i]    = static_cast<int>(d[i] - minfdelay);
            float fd = d[i] - static_cast<float>(id[i]);
            a[i]     = (1 - fd) / (1 + fd);
        }
        allpass_.setCoeff(a);
        TapOut::setDelay(id);
    }

    template <class Ctxt, class DL> auto read(Ctxt c, DL &delayline) const
    {
        auto x = TapOut::read(c, delayline);
        c.setIn(x);
        allpass_.process(c, delayline.getInner());
        return x;
    }

  private:
    AllPass<N, TapIn> allpass_;
};

template <int N> class TapAllPass<N, TapNoInterp<N>, TapTail> : TapNoInterp<N>
{
  public:
    static constexpr auto minfdelay = 0.618f;
    using TapOut                    = TapNoInterp<N>;
    using TapIn                     = TapTail;

    constexpr TapAllPass() = default;

    void setDelay(fData<N> d)
    {
        fData<N> a;
        iData<N> id;
        for (int i = 0; i < N; ++i) {
            id[i] = static_cast<int>(d[i] - minfdelay);

            float fd = d[i] - static_cast<float>(id[i]);
            // taylor approximation of (1-fd)/(1+fd)
            // -(fd-1)/2 + (fd-1)²/4 - (fd-1)³/8
            float z = (fd - 1.f) * 0.5f;
            a[i]    = z * (-1.f + z * (1.f - z));
        }
        allpass_.setCoeff(a);
        TapOut::setDelay(id);
    }

    template <class Ctxt, class DL> auto read(Ctxt c, DL &delayline) const
    {
        static_assert(Ctxt::VecSize == 1);
        typename Ctxt::Type x;
        typename Ctxt::Type x1;
        typename Ctxt::Type y;
        typename Ctxt::Type y1;

        for (int i = 0; i < N; i++) {
            auto id = TapOut::id_[i];
            assert(id + 1 <= DL::Length - Ctxt::VecSize + 1);

            auto &val = delayline.read(c, id);
            for (int k = 0; k < Ctxt::VecSize; ++k) x[k][i] = val[k][i];
            auto &val1 = delayline.read(c, id + 1);
            for (int k = 0; k < Ctxt::VecSize; ++k) x1[k][i] = val1[k][i];
        }

        y1 = allpass_.getTap().read(c, delayline.getInner());

        const auto &a = allpass_.getCoeff();
        for (int k = 0; k < Ctxt::VecSize; ++k) {
            for (int i = 0; i < N; i++) {
                y[k][i] = a[i] * (x[k][i] - y1[k][i]) + x1[k][i];
            }
        }
        delayline.getInner().write(c, y);
        return y;
    }

  private:
    AllPass<N, TapIn> allpass_;
};

template <int N> class AllPass2
{
    /* 2nd order allpass */

  public:
    class DL : public CopyDelayLine<N, 2>
    {
    };

    /**
     * @brief Set the poles of the allpass filter.
     *
     * This function sets the poles of the allpass filter based on the given
     * radius and frequency.
     *
     * @param R A signal representing the radius of the pole in the z-plane.
     *          Must be between -1 and 1.
     * @param freq A signal representing the frequency of the pole in the
     * z-plane. Typically in the range [0, 1], where 1 represents the Nyquist
     * frequency.
     */
    void setPole(fData<N> R, fData<N> freq)
    {
        auto &a1 = a_[1];
        auto &a2 = a_[0];
        for (int i = 0; i < N; ++i) {
            assert(R[i] >= -1 && R[i] <= 1);
            a2[i] = R[i] * R[i];
            a1[i] = -2 * R[i] * cos(M_PIf * freq[i]);
        }
    }

    template <class Ctxt, class DL> void process(Ctxt c, DL &delayline) const
    {
        // const auto &a1 = a_[1];
        const auto &a2 = a_[0];

        static_assert(Ctxt::VecSize == 1, "Vector size must be 1");

        // Input signal
        auto &x = c.getSignal();

        // Direct form II value and its delayed signals
        typename Ctxt::Type s0;
        auto &sN = delayline.read(c, 2)[0].toVector();
        // auto &s1 = sN[1];
        auto &s2 = sN[0];

        /* intermadiate values */
        typename Ctxt::BaseType sNa[2];

        // Apply the filter coefficients to the delayed signals.
        // Doing in a loop like this allows for compiler vectorization.
#pragma omp simd
        for (size_t k = 0; k < 2; ++k) {
            for (size_t i = 0; i < x[0].size(); ++i) {
                sNa[k][i] = sN[k][i] * a_[k][i];
            }
        }

        // Compute output using direct form II
#pragma omp simd
        for (size_t i = 0; i < x[0].size(); ++i) {
            s0[0][i] = x[0][i] - sNa[1][i] - sNa[0][i];
            x[0][i]  = a2[i] * s0[0][i] + sNa[1][i] + s2[i];
        }

        // Update the delay line with the new value
        delayline.write(c, s0);
    }

  private:
    fData<N> a_[2];
};
} // namespace dsp
