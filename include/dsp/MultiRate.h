#pragma once

#include "Buffer.h"
#include "FIRFilter.h"

namespace dsp
{

template <size_t N, size_t Order, size_t MaxM, size_t Offset = 0,
          size_t BufSize = 0, size_t M = MaxM>
class BaseMultiRate : public BaseMultiRate<N, Order, MaxM, Offset, BufSize, 1>
{
  public:
    using Base          = BaseMultiRate<N, Order, MaxM, Offset, BufSize, 1>;
    using Next          = BaseMultiRate<N, Order, MaxM, Offset, BufSize, M - 1>;
    using BufCtxt       = typename Next::BufCtxt;
    using Ctxt          = typename Next::Ctxt;
    using DLDecimate    = typename Next::DLDecimate;
    using DLInterpolate = typename Next::DLInterpolate;

    template <size_t BufSize_>
    using WithBuffer = BaseMultiRate<N, Order, MaxM, Offset, BufSize_, M>;

    BaseMultiRate(float cutoff = 1.f) :
        firdecimate_(cutoff), firinterpolate_(cutoff)
    {
    }

    int decimate(BufCtxt cin, Ctxt &cout, DLDecimate &dl,
                 int decimateId) const override
    {
        auto id = firdecimate_.decimate(cin.vec(), cout, dl, decimateId);
        return id;
    }
    int interpolate(Ctxt cin, Ctxt &cout, DLInterpolate &dl,
                    int interpolateId) const override
    {
        auto id = firinterpolate_.interpolate(cin, cout, dl, interpolateId);
        return id;
    }

  private:
    const FIRDecimate<N, Order, M> firdecimate_;
    const FIRInterpolate<N, Order, M> firinterpolate_;
};

template <size_t N, size_t Order, size_t MaxM, size_t Offset, size_t BufSize>
class BaseMultiRate<N, Order, MaxM, Offset, BufSize, 1>
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
    MultiRate(float cutoff = 1.f) : Next(cutoff), multirate_(cutoff) {}

    using Next = MultiRate<N, Order, MaxM, Offset, BufSize, M - 1>;
    using Base = typename Next::Base;
    template <size_t BufSize_>
    using WithBuffer = MultiRate<N, Order, MaxM, Offset, BufSize_, M>;
    [[nodiscard]] const typename Next::Base *get(int i) const
    {
        if (i >= static_cast<int>(M)) {
            return &multirate_;
        } else {
            return Next::get(i);
        }
    }

  private:
    const BaseMultiRate<N, Order, MaxM, Offset, BufSize, M> multirate_;
};

template <size_t N, size_t Order, size_t MaxM, size_t Offset, size_t BufSize>
class MultiRate<N, Order, MaxM, Offset, BufSize, 1>
{
  public:
    MultiRate(float cutoff = 1.f) { (void)cutoff; }

    using Base          = BaseMultiRate<N, Order, MaxM, Offset, BufSize, 1>;
    using BufCtxt       = typename Base::BufCtxt;
    using Ctxt          = typename Base::Ctxt;
    using DLDecimate    = typename Base::DLDecimate;
    using DLInterpolate = typename Base::DLInterpolate;
    [[nodiscard]] const Base *get(int i) const
    {
        (void)i;
        return &multirate_;
    }

  private:
    const Base multirate_;
};

} // namespace dsp
