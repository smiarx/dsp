#pragma once

#include "Delay.h"
#include "simd/multi.h"

namespace dsp
{
template <typename T> class PitchShift : TapCubic<T>
{
    using IntT = intType<T>;

  public:
    void setFreq(T freq)
    {
        assert(dsp::all(freq > 0 && freq <= 1));
        period_ = T(2) / freq;
    }
    void setShift(T shift)
    {
        assert(dsp::all(shift >= T(0)));
        shift_ = shift;
    }
    void setDelay(T delay)
    {
        assert(dsp::all(delay > T(1)));
        delay_ = delay;
    }

    template <class Ctxt, class DL> void process(Ctxt ctxt, DL delayline)
    {
        delayline.write(ctxt, ctxt.getInput());

        delay_ += T(1) - dsp::load(shift_);

        delay_ = dsp::blend(delay_ < T(2), delay_ + period_, delay_);
        delay_ = dsp::blend(delay_ > period_ + T(2), delay_ - period_, delay_);

        TapCubic<T>::setDelay(delay_);
        auto x = TapCubic<T>::read(ctxt, delayline);
        ctxt.setOutput(x);
    }

  private:
    T delay_{2};
    T period_{};
    T shift_{T(1)};
};
} // namespace dsp
