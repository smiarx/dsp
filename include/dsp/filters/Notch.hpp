#pragma once

#include "dsp/math/FastMath.hpp"
#include "dsp/simd/simd_math.hpp"
#include <array>
#include <dsp/cpu/defines.h>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T> class NotchFilter
{
  public:
    class State
    {
      public:
        friend NotchFilter;
        State(T a = T(1.99)) : a_{a} {}
        void setState(T value) { s_[0] = s_[1] = value; }

        const auto &getState() { return s_; }

        T getCoeff() const { return a_; }

        void setCoeff(T a) { a_ = a; }

        T getFreq() const
        {
            // get frequency from notch coeff
            return std::acos(a_ / 2) / dsp::constants<T>::pi;
        }

        void setFreq(T freq)
        {
            a_ = T(2) * dsp::cos(dsp::constants<T>::pi * freq);
        }

      private:
        T a_{1.999};           // notch parameter
        std::array<T, 2> s_{}; // state value
    };

    NotchFilter(T rho = T(0.95)) : rho_(rho), rhosq_(rho * rho) {}

    T getRho() const { return rho_; }

    void setRho(T rho)
    {
        rho_   = rho;
        rhosq_ = rho * rho;
    }

    template <class Ctxt> void process(Ctxt ctxt, State &state) const
    {
        auto x = ctxt.getInput();

        auto &a = state.a_;
        auto &s = state.s_;

        //  pole filter
        auto sn = x + rho_ * a * s[0] - rhosq_ * s[1];
        // zero filter
        auto y = sn - a * s[0] + s[1];

        ctxt.setOutput(y);
        s[1] = s[0];
        s[0] = sn;
    }

  private:
    T rho_{0.95}; // pole radius
    T rhosq_{rho_ * rho_};
};

/////////////////// Adaptive Notch Filter //////////////////////////
template <typename T> class AdaptiveNotchFilter : public NotchFilter<T>
{
    /* adaptive notch filter with recursive least square
     *
     * adapted from
     * https://www.dafx.de/paper-archive/2023/DAFx23_paper_7.pdf
     *
     * remove single pure frequency from signal, can also detect pitch
     * as notch cut frequency
     */

  public:
    class State : public NotchFilter<T>::State
    {
      public:
        friend AdaptiveNotchFilter;

        State(T a = T(1.99), T p = T(1)) : NotchFilter<T>::State(a), p_(p) {}

        void setFreq(T freq) { NotchFilter<T>::State::setFreq(freq); }

        void setFreq(T freq, T p)
        {
            setFreq(freq);
            p_ = p;
        }

        T getVariance() const { return p_; }

      private:
        T p_{1}; // variance
    };

    AdaptiveNotchFilter(T rho = T(0.95), T lambda = T(0.99)) :
        NotchFilter<T>(rho), lambda_(lambda), invlambda_(T(1) / lambda)
    {
    }

    T getLambda() const { return lambda_; }
    void setLambda(T lambda)
    {
        lambda_    = lambda;
        invlambda_ = T(1) / lambda;
    }

    template <class Ctxt> void process(Ctxt ctxt, State &state) const
    {
        NotchFilter<T>::process(ctxt, state);

        // error is notch output
        auto e = ctxt.getInput();

        auto a  = state.getCoeff();
        auto &s = state.getState();
        auto &p = state.p_;

        // kalman gain
        auto ps = p * s[1];
        auto k  = ps / (s[1] * ps + lambda_);

        // adjust
        a += k * e;
        p = invlambda_ * (p - k * ps);

        if (dsp::abs(a) > 2.0) {
            a = T(1.999);
        }

        state.setCoeff(a);
    }

  private:
    T lambda_{0.99}; // forgetting factor
    T invlambda_{1. / lambda_};
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
