#pragma once

#include "Delay.h"
#include "LinearAlgebra.h"
#include "dsp/AllPass.h"

namespace dsp
{

template <typename T, size_t Order> class AdaptiveFilter
{
  public:
    static constexpr auto kDelaySize = Order;
    using AnalyzeState               = CopyDelayLine<T, Order>;
    using ReconstructState           = CopyDelayLine<T, Order>;

    AdaptiveFilter() = default;
    AdaptiveFilter(const linalg::Vector<T, Order> &a) : a_(a) {}

    template <class Ctxt, class Delay>
    linalg::Vector<T, Order> analyze(Ctxt ctxt, Delay &delay) const
    {
        static_assert(!Ctxt::kUseVec);

        // x is previous inputs
        auto x = delay.asVector(ctxt);

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
        linalg::Vector<T, Order> y = delay.asVector(ctxt);
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

/////////////////////////////// Warped IIR Filter ////////////////////////

template <class T, size_t Order>
class WarpedAdaptiveFilter : public AdaptiveFilter<T, Order>
{
    // Implementation of frequency-warped recursive filters - Aki Härmä

    static constexpr auto kLength =
        nextAlignedOffset(Order, kUsedSIMDSize<T, Order>);

  public:
    // analyze delay line for FIR
    struct AnalyzeState {
        AllPassDelayLine<T, Order> delay_{};
        T compensateFilterMem_{};
    };
    // state for reconstruct IIR
    struct ReconstructState {
        std::array<T, Order> stateIIR{};
        T compensateFilterMem_{};
    };

    WarpedAdaptiveFilter() = default;
    WarpedAdaptiveFilter(const linalg::Vector<T, Order> &a) :
        AdaptiveFilter<T, Order>(a)
    {
    }

    void setWarpIn(T warpIn)
    {
        warpIn_     = warpIn;
        warpInGain_ = T(1) / std::sqrt(1 - warpIn * warpIn);
    }

    void setWarpOut(T warpOut)
    {
        warpOut_     = warpOut;
        warpOutGain_ = std::sqrt(1 - warpOut * warpOut);
    }

    template <class Ctxt, class State>
    void compensateResidualAnalyze(Ctxt ctxt, State &analyzeState) const
    {
        auto e = ctxt.getInput();

        // filter compensation of residual
        auto memE = e;
        e += warpIn_ * analyzeState.compensateFilterMem_;
        e *= warpInGain_;
        analyzeState.compensateFilterMem_ = memE;

        ctxt.setOutput(e);
    }

    template <class Ctxt, class State>
    void compensateResidualReconstruct(Ctxt ctxt, State &reconstructState) const
    {
        auto e = ctxt.getInput();
        // filter input with compensation filter
        e = e * warpOutGain_ - warpOut_ * reconstructState.compensateFilterMem_;
        reconstructState.compensateFilterMem_ = e;

        ctxt.setOutput(e);
    }

    template <class Ctxt, class State>
    linalg::Vector<T, Order> analyze(Ctxt ctxt, State &analyzeState) const
    {
        static_assert(!Ctxt::kUseVec);

        auto &a = AdaptiveFilter<T, Order>::getCoeffs();

        // set allpass delayline with warp coeff
        auto x = analyzeState.delay_.asVector(ctxt, warpIn_);

        // d is desired signal (current input)
        auto d = ctxt.getInput();

        // error
        auto e = d - a.dot(x);

        // write error in ouput
        ctxt.setOutput(e);

        // write new value to input
        analyzeState.delay_.write(ctxt, d);

        return x;
    }

    template <class Ctxt, class State>
    void reconstruct(Ctxt ctxt, State &reconstructState)
    {
        auto &a     = AdaptiveFilter<T, Order>::getCoeffs();
        auto &state = reconstructState.stateIIR;

        auto x = ctxt.getInput();

        // compute state and gain
        T in = state[0] - state[0] * warpOut_;
        T S  = a.get(Order - 1) * in;
        for (size_t n = 1; n < Order; ++n) {
            in = state[n] + warpOut_ * (in - state[n]);
            S += a.get(Order - 1 - n) * in;
        }

        T g = a.get(0) * (warpOut_);
        for (size_t n = 1; n < Order; ++n) {
            g += a.get(n);
            g *= warpOut_;
        }

        // output
        auto y = (x + S) / (1 - g);
        ctxt.setOutput(y);

        // update state
        auto out = y;
        for (size_t n = 0; n < Order; ++n) {
            auto in  = out;
            auto v   = warpOut_ * (in - state[n]);
            out      = v + state[n];
            state[n] = v + in;
        }
    }

  private:
    T warpIn_{};
    T warpInGain_{};
    T warpOut_{};
    T warpOutGain_{};
};

//////////////////////////// RLS //////////////

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

    template <class Ctxt, class AFilter, class State>
    void process(Ctxt ctxt, State &state, AFilter &filter)
    {
        // analyze with filter and write error
        auto x = filter.analyze(ctxt, state);
        auto e = ctxt.getInput();

        // update filter
        linalg::Vector<T, Order> phi, g;
        phi = P_.mul(x);
        g   = phi / (lambda_ + phi.dot(x));

        linalg::Matrix<T, Order, Order> out, outP;
        out  = g.outer(x);
        outP = out.mul(P_);
        P_   = (P_ - outP) * invLambda_;
        // we need do decompose this, why ?
        // P_     = (P_ - g.outer(x).mul(P_)) / lambda_;

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
    template <class Ctxt, class AFilter, class Delay>
    void process(Ctxt ctxt, Delay &delay, AFilter &filter)
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
