#pragma once

#include "Delay.h"
#include "LinearAlgebra.h"

namespace dsp
{

template <typename T, size_t Order> class AdaptiveFilter
{
  public:
    static constexpr auto kDelaySize = Order;
    using DL                         = CopyDelayLine<T, Order>;

    template <class Ctxt>
    void analyze(Ctxt ctxt, const linalg::Vector<T, Order> &x) const
    {
        static_assert(!Ctxt::kUseVec);

        // x is previous inputs
        // d is desired signal (current input)
        auto d = ctxt.getInput();

        // error
        auto e = d - a_.dot(x);

        // write error in ouput
        ctxt.setOutput(e);
    }

    void update(const linalg::Vector<T, Order> &deltaA) { a_ = a_ + deltaA; }

    const auto &getCoeffs() const { return a_; }

  private:
    linalg::Vector<T, Order> a_{}; // filter coefficients
};

template <typename T, size_t Order> class RLS
{
    // recursive least squares
  public:
    RLS(T lambda = 0.99)
    {
        setForgetFactor(lambda);

        for (size_t i = 0; i < Order; ++i) P_.set(i, i, 100);
    }

    void setForgetFactor(T lambda)
    {
        lambda_    = lambda;
        invLambda_ = 1 / lambda;
    }

    template <class Ctxt, class Delay>
    void process(Ctxt ctxt, Delay &delay, AdaptiveFilter<T, Order> &filter)
    {
        // get previous input as vector
        auto x = delay.asVector();

        // analyze with filter and write error
        filter.analyze(ctxt, x);
        auto e = ctxt.getInput();

        // update filter
        linalg::Vector<T, Order> phi, g;
        phi    = P_.mul(x);
        auto a = T{1} / (lambda_ + phi.dot(x));
        g      = phi * a;
        P_     = (P_ - g.outer(x).mul(P_)) * invLambda_;

        linalg::Vector<T, Order> deltaA;
        deltaA = g * e;
        filter.update(deltaA);
    }

  private:
    T lambda_{0.99};
    T invLambda_{1 / lambda_};
    linalg::Matrix<T, Order, Order> P_{}; // invert covariance matrix;
};

}; // namespace dsp
