#pragma once

#include "Buffer.h"
#include "FIRFilter.h"

namespace dsp
{

template <int N, int Order, int MaxM, int Offset = 0, int BufSize = 0,
          int M = MaxM>
class MultiRate : public MultiRate<N, Order, MaxM, Offset, BufSize, 1>
{
  public:
    using Base          = MultiRate<N, Order, MaxM, Offset, BufSize, 1>;
    using Next          = MultiRate<N, Order, MaxM, Offset, BufSize, M - 1>;
    using BufCtxt       = typename Next::BufCtxt;
    using Ctxt          = typename Next::Ctxt;
    using DLDecimate    = typename Next::DLDecimate;
    using DLInterpolate = typename Next::DLInterpolate;

    template <int BufSize_>
    using WithBuffer = MultiRate<N, Order, MaxM, Offset, BufSize_, M>;

    virtual void decimate(BufCtxt cin, Ctxt &cout, DLDecimate &dl)
    {
        firdecimate.decimate(cin.vec(), cout, dl);
    }
    virtual void interpolate(Ctxt cin, Ctxt &cout, DLInterpolate &dl)
    {
        firinterpolate.interpolate(cin, cout, dl);
    }

  private:
    FIRDecimate<N, Order, M> firdecimate;
    FIRInterpolate<N, Order, M> firinterpolate;
};

template <int N, int Order, int MaxM, int Offset, int BufSize>
class MultiRate<N, Order, MaxM, Offset, BufSize, 1>
{
  public:
    using BufCtxt = BufferContext<Signal<N>, Buffer<Signal<N>, BufSize>>;
    using Ctxt    = Context<Signal<N>>;
    using DLDecimate =
        typename FIRDecimate<N, Order, MaxM>::template DL<Offset>;
    using DLInterpolate =
        typename FIRInterpolate<N, Order, MaxM>::template DL<N>;

    virtual void decimate(BufCtxt cin, Ctxt &cout, DLDecimate &dl)
    {
        cout = cin;
    }
    virtual void interpolate(Ctxt cin, Ctxt &cout, DLInterpolate &dl)
    {
        cout = cin;
    }
};

template <int N, int Order, int MaxM, int Offset = 0, int BufSize = 0,
          int M = MaxM>
class MultiRates : public MultiRates<N, Order, MaxM, Offset, BufSize, M - 1>
{
  public:
    using Base = MultiRates<N, Order, MaxM, Offset, BufSize, 1>;
    using Next = MultiRates<N, Order, MaxM, Offset, BufSize, M - 1>;
    template <int BufSize_>
    using WithBuffer = MultiRates<N, Order, MaxM, Offset, BufSize_, M>;
    void set(unsigned int i)
    {
        if (i >= M) {
            Base::selected = &multirate;
        } else {
            Next::set(i);
        }
    }

  private:
    MultiRate<N, Order, MaxM, Offset, BufSize, M> multirate;
};

template <int N, int Order, int MaxM, int Offset, int BufSize>
class MultiRates<N, Order, MaxM, Offset, BufSize, 1>
{
  public:
    using Base          = MultiRate<N, Order, MaxM, Offset, BufSize, 1>;
    using BufCtxt       = typename Base::BufCtxt;
    using Ctxt          = typename Base::Ctxt;
    using DLDecimate    = typename Base::DLDecimate;
    using DLInterpolate = typename Base::DLInterpolate;
    void set(unsigned int i) { selected = &multirate; }

    virtual void decimate(BufCtxt cin, Ctxt &cout, DLDecimate &dl)
    {
        selected->decimate(cin, cout, dl);
    }
    virtual void interpolate(Ctxt cin, Ctxt &cout, DLInterpolate &dl)
    {
        selected->interpolate(cin, cout, dl);
    }

  protected:
    Base multirate;
    Base *selected{&multirate};
};

} // namespace dsp
