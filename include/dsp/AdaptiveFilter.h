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

    template <class Ctxt, class Delay>
    linalg::Vector<T, Order> analyze(Ctxt ctxt, Delay &delay) const
    {
        static_assert(!Ctxt::kUseVec);

        // x is previous inputs
        auto x = delay.asVector();

        // d is desired signal (current input)
        auto d = ctxt.getInput();

        // error
        auto e = d - a_.dot(x);

        // write error in ouput
        ctxt.setOutput(e);

        // write new value to input
        delay.write(ctxt, d);

        return x;
    }

    template <class Ctxt, class DL> void reconstruct(Ctxt ctxt, DL &delay)
    {
        linalg::Vector<T, Order> y = delay.asVector();
        auto e                     = ctxt.getInput();
        auto out                   = e + y.dot(a_);

        delay.write(ctxt, out);
        ctxt.setOutput(out);
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
        // analyze with filter and write error
        auto x = filter.analyze(ctxt, delay);
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

template <typename T, size_t Order> class RLSDCD
{
  public:
    RLSDCD(T lambda = 0.99)
    {
        setForgetFactor(lambda);

        for (size_t i = 0; i < Order; ++i) R_.set(i, i, 0.05);
    }

    void setForgetFactor(T lambda) { lambda_ = lambda; }

    // https://core.ac.uk/download/pdf/1145733.pdf
    template <class Ctxt, class Delay>
    void process(Ctxt ctxt, Delay &delay, AdaptiveFilter<T, Order> &filter)
    {
        // analyze with filter and write error
        auto x = filter.analyze(ctxt, delay);
        auto e = ctxt.getInput();

        // covariance matrix
        R_ = R_ * lambda_ + x.outer(x);

        // residual
        r_                = r_ * lambda_ + x * e;
        const auto deltaA = dcd();

        filter.update(deltaA);
    }

    auto dcd()
    {
        // dichotomous coordinate descent
        linalg::Vector<T, Order> deltaA{};
        T alpha = 2.f;
        int m   = 1;

        for (int k = 1; k < 128; ++k) {
            // arg max r
            int n   = 0;
            auto rn = r_.get(0);
            for (size_t i = 1; i < Order; ++i) {
                auto ri = r_.get(i);
                if (std::abs(ri) > std::abs(rn)) {
                    n  = i;
                    rn = ri;
                }
            }
            auto Rnn = R_.get(n, n);

            while (std::abs(rn) <= alpha / 2 * Rnn) {
                ++m;
                if (m == 28) {
                    return deltaA;
                }
                alpha = alpha / 2;
            }

            if (rn > 0) {
                deltaA.set(n, 0, deltaA.get(n) + alpha);
                r_ = r_ - R_.column(n) * alpha;
            } else {
                deltaA.set(n, 0, deltaA.get(n) - alpha);
                r_ = r_ + R_.column(n) * alpha;
            }
        }
        return deltaA;
    }

  private:
    T lambda_{0.99};
    linalg::Vector<T, Order> r_{};        // residual solution
    linalg::Matrix<T, Order, Order> R_{}; // covariance matrix;
};
}; // namespace dsp
