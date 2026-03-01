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
        kFadeIn,
        kOn,
        kFadeOut,
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

    void setPlayRate(T rate)
    {
        assert(rate > T(0));
        playRate_ = rate;
    }

    const auto &getFade() const { return fade_; }

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

        phase = -T(1) + synthPhase_ + marker.period_ * synthRate_ / playRate_;

        //  center position of grain
        auto grainPos = ctxt.getId() - phase * period;
        delay         = ctxt.modDiff(grainPos, marker.position_);

        if (delay < 0.f) return;

        // std::cout << 48000.f*rate << "\n";
        Grain newGrain{phase, rate, delay};

        grains_[nGrains_] = newGrain;
        nGrains_++;
    }

    template <class Ctxt, class Analysis>
    void start(Ctxt ctxt, Analysis analysis)
    {
        state_       = kFadeIn;
        fade_.value_ = std::max(-1.f, -fade_.value_);

        const auto &marker = analysis.getClosestMarker(ctxt);
        auto period        = marker.period_;
        // synthPhase_ = T(1) + ctxt.idDiff(marker.position_ + period) / period;
        synthPhase_ =
            T(1) + ctxt.idDiff(marker.position_ + period) /
                       period; // - marker.period_*synthRate_/playRate_;

        nGrains_ = 0;
        addGrain(ctxt, analysis);
        fade_.rate_ = grains_[0].rate_ * T(0.5);
    }

    template <class Ctxt> void stop(Ctxt ctxt)
    {
        state_       = kFadeOut;
        fade_.value_ = std::max(0.f, -fade_.value_);
        if (nGrains_ > 0) fade_.rate_ = grains_[0].rate_ * T(0.5);
        else {
            constexpr auto kDefaultFadeRate = T(0.02);
            fade_.rate_                     = kDefaultFadeRate;
        }
    }

    template <class Ctxt, class DL, class Analysis>
    void process(Ctxt ctxt, DL delayline, Analysis analysis)
    {
        if (state_ == kOff) {
            return;
        }

        const auto &marker = analysis.getClosestMarker(ctxt, playRate_);
        auto currentPhase =
            synthPhase_ + marker.period_ * synthRate_ / playRate_;
        if (currentPhase > T(1)) {
            synthPhase_ -= T(1);
            addGrain(ctxt, analysis);
        }

        auto x = T(0);
        for (size_t gId = 0; gId < nGrains_; ++gId) {
            auto &grain = grains_[gId];

            if (grain.phase_ > T(1)) {
                // delete grain (swap last)
                grains_[gId] = grains_[--nGrains_];
                --gId;
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

        T fadeValue;
        if (state_ == kFadeIn || state_ == kFadeOut) {
            auto &value = fade_.value_;
            value += fade_.rate_;
            if (state_ == kFadeIn && value > T(0)) {
                value  = T(0);
                state_ = kOn;
            } else if (state_ == kFadeOut && value > T(1)) {
                value  = T(1);
                state_ = kOff;
            }
            // TODO
            fadeValue = T(1) - Window::generate(value);
        } else {
            fadeValue = T(0);
        }

        x += fadeValue * (ctxt.getInput() - x);

        ctxt.setOutput(x);

        // increment synth phase
        synthPhase_ += synthRate_;
    }

  private:
    T synthRate_{1.f};

    // generating
    State state_{kOff};
    struct Fade {
        float value_;
        float rate_;
    } fade_;
    size_t nGrains_{};
    std::array<Grain, MaxGrains> grains_{};

    T synthPhase_{};
    T playRate_{1};
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
