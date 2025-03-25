#pragma once

#include "Buffer.h"
#include "FIRFilter.h"

namespace dsp
{

template <size_t N, size_t Order, size_t MaxM, size_t Offset = 0,
          size_t BufSize = 0, size_t M = MaxM>
class _MultiRate : public _MultiRate<N, Order, MaxM, Offset, BufSize, 1>
{
  public:
    using Base          = _MultiRate<N, Order, MaxM, Offset, BufSize, 1>;
    using Next          = _MultiRate<N, Order, MaxM, Offset, BufSize, M - 1>;
    using BufCtxt       = typename Next::BufCtxt;
    using Ctxt          = typename Next::Ctxt;
    using DLDecimate    = typename Next::DLDecimate;
    using DLInterpolate = typename Next::DLInterpolate;

    template <size_t BufSize_>
    using WithBuffer = _MultiRate<N, Order, MaxM, Offset, BufSize_, M>;

    _MultiRate(float cutoff = 1.f) : firdecimate(cutoff), firinterpolate(cutoff)
    {
    }

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

template <size_t N, size_t Order, size_t MaxM, size_t Offset, size_t BufSize>
class _MultiRate<N, Order, MaxM, Offset, BufSize, 1>
{
  public:
    using BufCtxt = BufferContext<fSample<N>, Buffer<fSample<N>, BufSize>>;
    using Ctxt    = Context<fSample<N>>;
    using DLDecimate =
        typename FIRDecimate<N, Order, MaxM>::template DL<Offset>;
    using DLInterpolate =
        typename FIRInterpolate<N, Order, MaxM>::template DL<N>;

    virtual int decimate(BufCtxt cin, Ctxt &cout, DLDecimate &dl, int) const
    {
        (void)dl;
        cout = cin;
        return 0;
    }
    virtual int interpolate(Ctxt cin, Ctxt &cout, DLInterpolate &dl, int) const
    {
        (void)dl;
        cout = cin;
        return 0;
    }
};

template <size_t N, size_t Order, size_t MaxM, size_t Offset = 0,
          size_t BufSize = 0, size_t M = MaxM>
class MultiRate : public MultiRate<N, Order, MaxM, Offset, BufSize, M - 1>
{
  public:
    MultiRate(float cutoff = 1.f) : Next(cutoff), multirate(cutoff) {}

    using Next = MultiRate<N, Order, MaxM, Offset, BufSize, M - 1>;
    using Base = typename Next::Base;
    template <size_t BufSize_>
    using WithBuffer = MultiRate<N, Order, MaxM, Offset, BufSize_, M>;
    const typename Next::Base *get(int i) const
    {
        if (i >= static_cast<int>(M)) {
            return &multirate;
        } else {
            return Next::get(i);
        }
    }

  private:
    const _MultiRate<N, Order, MaxM, Offset, BufSize, M> multirate;
};

template <size_t N, size_t Order, size_t MaxM, size_t Offset, size_t BufSize>
class MultiRate<N, Order, MaxM, Offset, BufSize, 1>
{
  public:
    MultiRate(float cutoff = 1.f) { (void)cutoff; }

    using Base          = _MultiRate<N, Order, MaxM, Offset, BufSize, 1>;
    using BufCtxt       = typename Base::BufCtxt;
    using Ctxt          = typename Base::Ctxt;
    using DLDecimate    = typename Base::DLDecimate;
    using DLInterpolate = typename Base::DLInterpolate;
    const Base *get(int i) const
    {
        (void)i;
        return &multirate;
    }

  private:
    const Base multirate;
};

} // namespace dsp
