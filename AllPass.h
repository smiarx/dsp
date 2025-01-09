#pragma once

#include "Delay.h"
#include "Signal.h"

namespace dsp
{
template <class... T> class List
{
};

template <int N, class Tap> class TapAllPass;

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
    friend TapAllPass<N, Tap>;

  public:
    constexpr AllPass(const Signal<N> &a) : a_(a) {}

    template <class Ctxt, class DL> void process(Ctxt c, DL &delayline) const
    {
        auto &x = c.getIn();

        auto sN = tap_.read(c, delayline);
        decltype(sN) s0;
        for (int k = 0; k < c.getStep(); ++k) {
            for (int i = 0; i < N; ++i) {
                s0[k][i] = x[k][i] - a_[i] * sN[k][i];
                x[k][i]  = a_[i] * s0[k][i] + sN[k][i];
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

template <int N, class Tap = TapFix<>> class TapAllPass : public TapNoInterp<N>
{
  public:
    static constexpr auto minfdelay = 0.1f;

    void setDelay(int i, float d)
    {
        auto id        = static_cast<int>(d - minfdelay);
        float fd       = d - static_cast<float>(id);
        allpass_.a_[i] = (1 - fd) / (1 + fd);
        TapNoInterp<N>::setDelay(i, id);
    };

    void setDelay(Signal<N> d)
    {
        for (int i = 0; i < N; ++i) {
            setDelay(i, d[i]);
        }
    }

    template <class Ctxt, class DL> Signal<N> read(Ctxt c, DL &delayline) const
    {
        auto x = TapNoInterp<N>::read(c, delayline);
        c.setIn(&x);
        allpass_.process(c, delayline.getInner());
        return x;
    }

  private:
    AllPass<N, Tap> allpass_;
};
} // namespace dsp
