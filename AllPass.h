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
    constexpr AllPass(const Signal<N> &a, Args... args) : a_(a), tap_(args...)
    {
    }

    template <class Ctxt, class DL> void process(Ctxt c, DL &delayline) const
    {
        static_assert(Ctxt::VecSize <= DL::Length);

        auto &x = c.getIn();
        auto sN = tap_.read(c, delayline);
        decltype(sN) s0;

        for (int k = 0; k < x.size(); ++k) {
            for (int i = 0; i < x[0].size(); ++i) {
                auto a = a_[i % N];

                s0[k][i] = x[k][i] - a * sN[k][i];
                x[k][i]  = a * s0[k][i] + sN[k][i];
            }
        }
        delayline.write(c, s0);
    }

    void setCoeff(Signal<N> a) { a_ = a; }
    void setDelay(Signal<N> d) { tap_.setDelay(d); }
    void setDelay(iSignal<N> d) { tap_.setDelay(d); }

  private:
    Signal<N> a_;
    Tap tap_;
};

template <int N, class TapOut = TapTail, class TapIn = TapTail>
class TapAllPass : public TapOut
{
  public:
    static constexpr auto minfdelay = 0.1f;

    constexpr TapAllPass() = default;
    template <class... Args>
    constexpr TapAllPass(const Signal<N> &a, Args... args) :
        TapOut(args...), allpass_(a)
    {
    }

    void setDelay(int i, float d)
    {
        auto id        = static_cast<int>(d - minfdelay);
        float fd       = d - static_cast<float>(id);
        allpass_.a_[i] = (1 - fd) / (1 + fd);
        TapOut::setDelay(i, id);
    };

    void setDelay(Signal<N> d)
    {
        for (int i = 0; i < N; ++i) {
            setDelay(i, d[i]);
        }
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

template <int N> class AllPass2
{
    /* 2nd order allpass */

  public:
    class DL : public CopyDelayLine<N, 2>
    {
    };

    void setPole(Signal<N> R, Signal<N> freq)
    {
        auto &a1 = a_[1];
        auto &a2 = a_[0];
        for (int i = 0; i < N; ++i) {
            a2[i] = R[i] * R[i];
            a1[i] = -2 * R[i] * cos(M_PIf * freq[i]);
        }
    }

    template <class Ctxt, class DL> void process(Ctxt c, DL &delayline) const
    {
        // const auto &a1 = a_[1];
        const auto &a2 = a_[0];

        static_assert(Ctxt::VecSize == 1);
        auto &x  = c.getIn();
        auto &sN = delayline.read(c, 2)[0].toVector();
        typename Ctxt::Type s0;
        // auto &s1 = sN[1];
        auto &s2 = sN[0];

        typename Ctxt::BaseType sNa[2];

#pragma omp simd
        for (size_t k = 0; k < 2; ++k) {
            for (size_t i = 0; i < x[0].size(); ++i) {
                sNa[k][i] = sN[k][i] * a_[k][i];
            }
        }

#pragma omp simd
        for (size_t i = 0; i < x[0].size(); ++i) {
            s0[0][i] = x[0][i] - sNa[1][i] - sNa[0][i];
            x[0][i]  = a2[i] * s0[0][i] + sNa[1][i] + s2[i];
        }
        delayline.write(c, s0);
    }

  private:
    Signal<N> a_[2];
};
} // namespace dsp
