#pragma once

#include <dsp/analysis/ZFF.hpp>
#include <dsp/core/Buffer.hpp>
#include <dsp/cpu/infos.h>
#include <dsp/delay/Delay.hpp>
#include <dsp/resampling/Windows.hpp>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T, size_t MaxGrains = 6, class Window = windows::Hann>
class PSOLA
{
  public:
    using position_t                 = dsp::bufferId_t;
    static constexpr auto kSearchMin = T(0.98);
    static constexpr auto kSearchMax = T(1.02);

    static constexpr int kMaxMarkers = 4;

    enum State {
        kOff,
        kWaitingForLastGrains,
        kOn,
    };

    struct Grain {
        T phase_{};
        T rate_{};
        T delay_{};
    };

    PSOLA() = default;

    void setFreq(T freq)
    {
        assert(freq > T(0));
        synthRate_ = freq * T(0.5);
    }

    void setShift(T shift) { shift_ = shift; }

    void setPlayRate(T rate)
    {
        assert(rate > T(0));
        playRate_ = rate;
    }

    auto getState() const { return state_; }

    auto getNumGrains() const { return nGrains_; }

    auto getGrain(size_t id) const
    {
        assert(id < static_cast<size_t>(nGrains_));
        return grains_[id];
    }

    template <class Ctxt, class Analysis>
    void addGrain(Ctxt ctxt, Analysis analysis)
    {
        if (nGrains_ >= MaxGrains) return;
        assert(nGrains_ < MaxGrains);

        auto delay = T(0);
        T phase, rate;

        const auto &marker = analysis.getClosestMarker(ctxt, playRate_);
        auto period        = marker.period_;

        assert(period > T(0));
        rate = T(1) / period;

        phase = -T(1) + synthPhase_ + shift_ / playRate_;

        //  center position of grain
        auto grainPos = static_cast<T>(ctxt.getId()) - phase * period;
        delay         = ctxt.modDiff(grainPos, marker.position_);

        if (delay < T(0)) return;

        Grain newGrain{phase, rate, delay};

        grains_[nGrains_] = newGrain;
        nGrains_++;
    }

    template <class Ctxt, class Analysis>
    void start(Ctxt ctxt, Analysis analysis)
    {
        const auto &marker = analysis.getClosestMarker(ctxt);
        auto period        = marker.period_;

        if (state_ == kOff) {
            auto diffMarker = ctxt.idDiff(marker.position_);
            synthPhase_     = (diffMarker - period) * synthRate_ * shift_;

            // ensure phase is > -1
            while (synthPhase_ < T(-1)) synthPhase_ += T(1);
        }
        state_ = kOn;
    }

    void stop() { state_ = kWaitingForLastGrains; }

    template <class Ctxt, class DL, class Analysis>
    void process(Ctxt ctxt, DL delayline, Analysis analysis)
    {
        delayline.write(ctxt, ctxt.getInput());

        if (state_ == kOff) {
            ctxt.setOutput(T(0));
            return;
        }

        auto currentPhase = synthPhase_ + shift_ / playRate_;
        if (currentPhase > T(1)) {
            synthPhase_ -= T(1);
            if (state_ == kOn) {
                addGrain(ctxt, analysis);
            }
        }

        auto x = T(0);
        for (size_t gId = 0; gId < nGrains_; ++gId) {
            auto &grain = grains_[gId];

            if (grain.phase_ > T(1)) {
                // delete grain (swap last)
                grains_[gId] = grains_[--nGrains_];
                --gId;

                // stop if waiting for last grains
                if (nGrains_ == 0 && state_ == kWaitingForLastGrains) {
                    state_ = kOff;
                }

                continue;
            }

            // TODO
            auto fadeValue = Window::generate(grain.phase_);
            T grainValue;
            if (grain.delay_ > 2.f) {
                grainValue = TapCubic<T>{grain.delay_}.read(ctxt, delayline);
            } else {
                grainValue = TapLin<T>{grain.delay_}.read(ctxt, delayline);
            }

            x += grainValue * fadeValue;

            grain.phase_ += grain.rate_ * playRate_;
            grain.delay_ += T(1) - playRate_;
        }

        ctxt.setOutput(x);

        // increment synth phase
        synthPhase_ += synthRate_ * shift_;
    }

  private:
    T synthRate_{1.f};

    // generating
    State state_{kOff};
    size_t nGrains_{};
    std::array<Grain, MaxGrains> grains_{};

    T synthPhase_{};
    T playRate_{1};

    T shift_{1};
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
