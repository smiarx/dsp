#pragma once

#include <dsp/cpu/defines.h>
#include <dsp/simd/multi.hpp>
#include <dsp/simd/simd_math.hpp>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T> class Amplitude
{
  public:
    using State = T;

    Amplitude() = default;

    void setRate(T target, T up, T down)
    {
        upRate_   = dsp::pow(dsp::load(target), T(1) / up);
        downRate_ = dsp::pow(dsp::load(target), T(1) / down);
    }

    template <class Ctxt> void process(Ctxt ctxt, State &state)
    {
        auto in = ctxt.getInput();

        auto mask = in < state;
        auto rate = dsp::blend(mask, dsp::load(downRate_), dsp::load(upRate_));

        state = in + rate * (state - in);

        ctxt.setOutput(state);
    }

  private:
    T downRate_{}, upRate_{};
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
