#pragma once

#include <dsp/delay/Delay.hpp>
#include <dsp/resampling/Kernel.hpp>
#include <dsp/resampling/Windows.hpp>
#include <dsp/simd/multi.hpp>

namespace dsp
{
template <typename T> class PitchShift
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
        delayline.writeSafe(ctxt, ctxt.getInput());

        delay_ += T(1) - dsp::load(shift_);

        delay_ = dsp::blend(delay_ < T(5), delay_ + period_, delay_);
        delay_ = dsp::blend(delay_ > period_ + T(5), delay_ - period_, delay_);

        Tap tap;
        tap.setDelay(delay_);
        auto x = tap.read(ctxt, delayline, shift_);
        // auto x = tap.read(ctxt, delayline);
        ctxt.setOutput(x);
    }

  private:
    T delay_{2};
    T period_{};
    T shift_{T(1)};

    using Tap =
        dsp::TapKernel<T, dsp::kernels::Sinc<3, dsp::windows::Kaiser<140>>, 64>;
};
} // namespace dsp
