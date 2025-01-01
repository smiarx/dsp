#pragma once

#include "Delay.h"
#include "Signal.h"

namespace dsp
{
template <class... T> class List
{
};

template <int N, class Tap> class TapAllPass;

template <int N, class Tap = TapFix<N>> class AllPass
{
    friend TapAllPass<N, Tap>;

  public:
    template <int BufferSize, class D, class... Ds>
    void process(Context<N, BufferSize> c, Signal<N> &__restrict x,
                 D &delayline, Ds &...ds)
    {
        Signal<N> s0;
        Signal<N> sN = tap_.read(c, delayline, ds...);
        for (int i = 0; i < N; ++i) {
            s0[i] = x[i] - a_[i] * sN[i];
            x[i]  = a_[i] * s0[i] + sN[i];
        }
        delayline.write(c, s0);
    }

    void setCoeff(Signal<N> a) { a_ = a; }

    void setDelay(Signal<N> d) { tap_.setDelay(d); }

  private:
    Signal<N> a_;
    Tap tap_;
};

template <int N, class Tap = TapFix<N>> class TapAllPass : public TapNoInterp<N>
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

    template <int BufferSize, class D, class... Ds>
    Signal<N> read(Context<N, BufferSize> c, D &d, Ds &...ds)
    {
        auto x = TapNoInterp<N>::read(c, d);
        allpass_.process(c, x, ds...);
        return x;
    }

  private:
    AllPass<N, Tap> allpass_;
};
} // namespace dsp
