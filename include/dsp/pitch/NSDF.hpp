#pragma once

#include <cstddef>
#include <dsp/cpu/infos.h>
#include <dsp/delay/Delay.hpp>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T, size_t maxTau, size_t minTau = 2> class NSDF
{
    /* normalised square difference function */
  public:
    static constexpr auto kMinTau      = minTau;
    static constexpr auto kMaxTau      = maxTau;
    static constexpr size_t kTauLength = maxTau - minTau;

    static constexpr auto kMaxPeaks = kTauLength / 3;

    template <size_t Offset = 0> using DL = DelayLine<2 * maxTau, Offset>;

    void clear()
    {
        std::fill(std::begin(ac_), std::end(ac_), T(0));
        std::fill(std::begin(energy_), std::end(energy_), T(0));
    }

    T get(size_t tau) const
    {
        assert(tau >= kMinTau && tau < kMaxTau);

        auto i = tau2id(tau);
        return T(2) * ac_[i] / energy_[i];
    }

    auto updatePeaks(T threshold = T(0.01))
    {
        size_t nPeaks = 0;

        bool firstPeak            = true;
        bool negativeZeroCrossing = false;

        size_t localPeakId = 0;
        T peakAc           = T(0);
        T peakEnergy       = T(1000000);

        for (size_t tau = minTau; tau < maxTau; ++tau) {
            auto i  = tau2id(tau);
            auto ac = ac_[i];
            if (ac >= 0) {
                negativeZeroCrossing = false;
                if (!firstPeak) {
                    // auto value = ac_[i] / energy_[i];
                    if (ac_[i] * peakEnergy > peakAc * energy_[i]) {
                        localPeakId = i;
                        peakAc      = ac_[i];

                        peakEnergy = energy_[i];
                    }
                } else
                    continue;
            } else {
                if (!negativeZeroCrossing && localPeakId > 0) {
                    negativeZeroCrossing = true;
                    auto i               = localPeakId;
                    auto peakE = energy_[i] / static_cast<T>(id2tau(i));
                    if (peakE > threshold) {
                        auto newPeak =
                            energy_[i + 1] > T(0) && energy_[i - 1] > T(0)
                                ? extrema(interp::LagrangeParabolic{
                                      getBin(i + 1), getBin(i), getBin(i - 1)})
                                : interp::PolyPoint<T>{0, ac_[i] / energy_[i]};

                        assert(newPeak.x > -1.f && newPeak.x < 1.f);
                        newPeak.x += static_cast<T>(id2tau(i));

                        peaks_[nPeaks] = newPeak;
                        nPeaks++;
                        if (nPeaks >= kMaxPeaks) break;
                    }
                    peakAc     = T(0);
                    peakEnergy = T(1000000);
                }
                firstPeak = false;
                continue;
            }
        }

        return nPeaks;
    }

    [[nodiscard]] auto getBin(size_t id) { return ac_[id] / energy_[id]; }

    // should use a span be we still in c++17
    const auto &getPeaks() const { return peaks_; }

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
            auto xt  = delayline.read(ctxt, static_cast<int>(tau));
            auto x2t = delayline.read(ctxt, static_cast<int>(2 * tau));

            // if vectorized we need to get double tau vector
            // [t,t-1,t-2,t-3]*2 = [t,t-2,t-4,t-6]
            // we use the even function to get the values
            if constexpr (decltype(ctxt)::kIncrSize > 1) {
                auto x2t2 = delayline.read(
                    vecCtxt, static_cast<int>(2 * tau - kIncrSize));
                x2t = dsp::even(x2t, x2t2);
            }
            e += xsq - x2t * x2t;
            ac += (x - x2t) * xt;

            using eT = decltype(e);
            // prevent energy values from begin under 0;
            e = dsp::blend(e < eT(0), eT(0), e);
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

    interp::PolyPoint<T> peaks_[kMaxPeaks];

    T lastPeak_{};
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
