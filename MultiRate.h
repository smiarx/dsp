#pragma once

#include "Buffer.h"
#include "FIRFilter.h"

namespace dsp
{

template <int N, int Order, int MaxM, int Offset = 0, int BufSize = 0,
          int M = MaxM>
class _MultiRate : public _MultiRate<N, Order, MaxM, Offset, BufSize, 1>
{
  public:
    using Base          = _MultiRate<N, Order, MaxM, Offset, BufSize, 1>;
    using Next          = _MultiRate<N, Order, MaxM, Offset, BufSize, M - 1>;
    using BufCtxt       = typename Next::BufCtxt;
    using Ctxt          = typename Next::Ctxt;
    using DLDecimate    = typename Next::DLDecimate;
    using DLInterpolate = typename Next::DLInterpolate;

    template <int BufSize_>
    using WithBuffer = _MultiRate<N, Order, MaxM, Offset, BufSize_, M>;

    virtual int decimate(BufCtxt cin, Ctxt &cout, DLDecimate &dl,
                         int decimateId) const
    {
        auto id = firdecimate.decimate(cin.vec(), cout, dl, decimateId);
        return id;
    }
    virtual int interpolate(Ctxt cin, Ctxt &cout, DLInterpolate &dl,
                            int interpolateId) const
    {
        auto id = firinterpolate.interpolate(cin, cout, dl, interpolateId);
        return id;
    }

  private:
    const FIRDecimate<N, Order, M> firdecimate;
    const FIRInterpolate<N, Order, M> firinterpolate;
};

template <int N, int Order, int MaxM, int Offset, int BufSize>
class _MultiRate<N, Order, MaxM, Offset, BufSize, 1>
{
  public:
    using BufCtxt = BufferContext<Signal<N>, Buffer<Signal<N>, BufSize>>;
    using Ctxt    = Context<Signal<N>>;
    using DLDecimate =
        typename FIRDecimate<N, Order, MaxM>::template DL<Offset>;
    using DLInterpolate =
        typename FIRInterpolate<N, Order, MaxM>::template DL<N>;

    virtual int decimate(BufCtxt cin, Ctxt &cout, DLDecimate &dl, int) const
    {
        cout = cin;
        return 0;
    }
    virtual int interpolate(Ctxt cin, Ctxt &cout, DLInterpolate &dl, int) const
    {
        cout = cin;
        return 0;
    }
};

template <int N, int Order, int MaxM, int Offset = 0, int BufSize = 0,
          int M = MaxM>
class MultiRate : public MultiRate<N, Order, MaxM, Offset, BufSize, M - 1>
{
  public:
    MultiRate() = default;

    using Next = MultiRate<N, Order, MaxM, Offset, BufSize, M - 1>;
    using Base = typename Next::Base;
    template <int BufSize_>
    using WithBuffer = MultiRate<N, Order, MaxM, Offset, BufSize_, M>;
    const typename Next::Base *get(unsigned int i) const
    {
        if (i >= M) {
            return &multirate;
        } else {
            return Next::get(i);
        }
    }

  private:
    const _MultiRate<N, Order, MaxM, Offset, BufSize, M> multirate;
};

template <int N, int Order, int MaxM, int Offset, int BufSize>
class MultiRate<N, Order, MaxM, Offset, BufSize, 1>
{
  public:
    MultiRate() = default;

    using Base          = _MultiRate<N, Order, MaxM, Offset, BufSize, 1>;
    using BufCtxt       = typename Base::BufCtxt;
    using Ctxt          = typename Base::Ctxt;
    using DLDecimate    = typename Base::DLDecimate;
    using DLInterpolate = typename Base::DLInterpolate;
    const Base *get(unsigned int i) const { return &multirate; }

  private:
    const Base multirate;
};

} // namespace dsp
