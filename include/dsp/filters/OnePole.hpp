#pragma once

#include <dsp/core/Units.hpp>
#include <dsp/math/FastMath.hpp>
#include <dsp/math/math.hpp>

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

    void setConvergenceTime(T convergenceDiff, time::Samples<T> samples)
    {
        assert(dsp::all(convergenceDiff >= T(0) && convergenceDiff <= T(1) &&
                        samples.get() > T(0)));
        coeff_ = dsp::pow(dsp::load(convergenceDiff),
                          T(1) / dsp::load(samples.get()));
    }
    void setCoeff(T coeff)
    {
        assert(coeff >= T(0) && coeff <= T(1));
        coeff_ = coeff;
    }
    void setCutoff(freq::NyquistRatio<T> freq)
    {
        assert(freq >= T(0) && freq <= T(1));
        coeff_ = dsp::exp(-dsp::constants<T>::pi * dsp::load(freq.get()));
    }
    auto getCoeff() const { return coeff_; }

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
