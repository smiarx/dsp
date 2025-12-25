#pragma once

#include "Delay.h"
#include <cstddef>

namespace dsp
{

template <typename T, size_t maxTau, size_t minTau = 2> class NSDF
{
    /* normalised square difference function */
  public:
    static constexpr auto kMinTau      = minTau;
    static constexpr auto kMaxTau      = maxTau;
    static constexpr size_t kTauLength = maxTau - minTau;

    template <size_t Offset = 0> using DL = DelayLine<2 * maxTau, Offset>;

    T get(size_t tau) const
    {
        assert(tau >= kMinTau && tau < kMaxTau);

        auto i = tau2id(tau);
        return T(2) * ac_[i] / energy_[i];
    }

    template <class Ctxt, class DL> void process(Ctxt ctxt, DL &delayline)
    {
        auto x   = ctxt.getInput();
        auto xsq = x * x;

        auto vecCtxt             = ctxt.vec();
        constexpr auto kIncrSize = decltype(vecCtxt)::kIncrSize;

        auto nsdfId = [&](auto ctxt, size_t i) {
            auto tau = id2tau(i);
            auto e   = ctxt.load(energy_[i]);
            auto ac  = ctxt.load(ac_[i]);
            auto xt  = delayline.read(ctxt, tau);
            auto x2t = delayline.read(ctxt, 2 * tau);

            // if vectorized we need to get double tau vector
            // [t,t-1,t-2,t-3]*2 = [t,t-2,t-4,t-6]
            // we use the even function to get the values
            if constexpr (decltype(ctxt)::kIncrSize > 1) {
                auto x2t2 = delayline.read(vecCtxt, 2 * tau - kIncrSize);
                x2t       = dsp::even(x2t, x2t2);
            }
            e += xsq - x2t * x2t;
            ac += (x - x2t) * xt;

            ctxt.store(energy_[i], e);
            ctxt.store(ac_[i], ac);
        };

        size_t i;
        // do vectorized
        for (i = 0; i < kTauLength - kIncrSize + 1; i += kIncrSize) {
            nsdfId(vecCtxt, i);
        }
        // do scalar
        for (; i < kTauLength; ++i) {
            nsdfId(ctxt, i);
        }

        delayline.writeSafe(ctxt, x);
    }

  private:
    static constexpr auto tau2id(size_t tau) { return kMaxTau - 1 - tau; }
    static constexpr auto id2tau(size_t id) { return kMaxTau - 1 - id; }

    std::array<T, kTauLength> energy_{};
    std::array<T, kTauLength> ac_{};
};
} // namespace dsp
