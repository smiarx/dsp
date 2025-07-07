#pragma once

#include "Buffer.h"
#include "FIRFilter.h"

namespace dsp
{

template <typename T, size_t Order, size_t MaxM, size_t DecimOffset = 0,
          size_t InterpOffset = 0, size_t BufSize = 0, size_t M = MaxM>
class BaseMultiRate : public BaseMultiRate<T, Order, MaxM, DecimOffset,
                                           InterpOffset, BufSize, 1>
{
  public:
    using Base =
        BaseMultiRate<T, Order, MaxM, DecimOffset, InterpOffset, BufSize, 1>;
    using Next       = BaseMultiRate<T, Order, MaxM, DecimOffset, InterpOffset,
                                     BufSize, M - 1>;
    using BufCtxt    = typename Next::BufCtxt;
    using Ctxt       = typename Next::Ctxt;
    using DLDecimate = typename Next::DLDecimate;
    using DLInterpolate = typename Next::DLInterpolate;

    template <size_t BufSize_>
    using WithBuffer =
        BaseMultiRate<T, Order, MaxM, DecimOffset, InterpOffset, BufSize_, M>;

    BaseMultiRate(baseType<T> cutoff = 1) :
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
    const FIRDecimate<T, Order, M> firdecimate_;
    const FIRInterpolate<T, Order, M> firinterpolate_;
};

template <typename T, size_t Order, size_t MaxM, size_t DecimOffset,
          size_t InterpOffset, size_t BufSize>
class BaseMultiRate<T, Order, MaxM, DecimOffset, InterpOffset, BufSize, 1>
{
  public:
    using BufCtxt = BufferContext<T, BufSize>;
    using Ctxt    = Context<T>;
    using DLDecimate =
        typename FIRDecimate<T, Order, MaxM>::template DL<DecimOffset>;
    using DLInterpolate =
        typename FIRInterpolate<T, Order, MaxM>::template DL<InterpOffset>;

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

template <typename T, size_t Order, size_t MaxM, size_t DecimOffset = 0,
          size_t InterpOffset = 0, size_t BufSize = 0, size_t M = MaxM>
class MultiRate : public MultiRate<T, Order, MaxM, DecimOffset, InterpOffset,
                                   BufSize, M - 1>
{
  public:
    MultiRate(baseType<T> cutoff = 1) : Next(cutoff), multirate_(cutoff) {}

    using Next =
        MultiRate<T, Order, MaxM, DecimOffset, InterpOffset, BufSize, M - 1>;
    using Base = typename Next::Base;
    template <size_t BufSize_>
    using WithBuffer =
        MultiRate<T, Order, MaxM, DecimOffset, InterpOffset, BufSize_, M>;
    [[nodiscard]] const typename Next::Base *get(int i) const
    {
        if (i >= static_cast<int>(M)) {
            return &multirate_;
        } else {
            return Next::get(i);
        }
    }

  private:
    const BaseMultiRate<T, Order, MaxM, DecimOffset, InterpOffset, BufSize, M>
        multirate_;
};

template <typename T, size_t Order, size_t MaxM, size_t DecimOffset,
          size_t InterpOffset, size_t BufSize>
class MultiRate<T, Order, MaxM, DecimOffset, InterpOffset, BufSize, 1>
{
  public:
    MultiRate(float cutoff = 1.f) { (void)cutoff; }

    using Base =
        BaseMultiRate<T, Order, MaxM, DecimOffset, InterpOffset, BufSize, 1>;
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
