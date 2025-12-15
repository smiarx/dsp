#pragma once

#include "dsp/FastMath.h"
#include "dsp/simd/simd_math.h"

#include <cassert>

namespace dsp
{

template <typename T> class OnePole
{
    /*
     *   y[n] = (1-a)⋅x[n] + a⋅y[n-1]
     */
  public:
    using State = T;

    void setRate(T convergenceDiff, T rate)
    {
        assert(convergenceDiff >= T(0) && convergenceDiff <= T(1) &&
               rate >= T(0));
        coeff_ = dsp::pow(dsp::load(convergenceDiff), dsp::load(rate));
    }
    void setCoeff(T coeff)
    {
        assert(coeff >= T(0) && coeff <= T(1));
        coeff_ = coeff;
    }
    void setCutoff(T freq)
    {
        assert(freq >= T(0) && freq <= T(1));
        coeff_ = dsp::exp(-dsp::constants<T>::pi * dsp::load(freq));
    }
    auto getCoeff() const { return coeff_;}

    template <class Ctxt> void process(Ctxt ctxt, State &state) const
    {
        auto x = ctxt.getInput();

        x += coeff_ * (state - x);
        state = x;

        ctxt.setOutput(x);
    }

  private:
    T coeff_;
};
} // namespace dsp
