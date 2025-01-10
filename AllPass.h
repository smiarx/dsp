#pragma once

#include "Delay.h"
#include "Signal.h"

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
} // namespace dsp
