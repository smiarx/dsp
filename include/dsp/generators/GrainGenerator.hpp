#pragma once

#include "dsp/delay/Delay.hpp"
#include <cstddef>
#include <dsp/core/Buffer.hpp>
#include <dsp/cpu/defines.h>
#include <dsp/resampling/Windows.hpp>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T, size_t MaxGrains = 6, class Window = windows::Hann>
class GrainGenerator
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

    GrainGenerator() = default;

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

    auto getState() const { return state_; }

    auto getNumGrains() const { return nGrains_; }

    auto getGrain(size_t id) const
    {
        assert(id < static_cast<size_t>(nGrains_));
        return grains_[id];
    }

    void addGrain(Grain &&grain)
    {
        if (nGrains_ >= MaxGrains) return;
        grains_[nGrains_] = grain;
        nGrains_++;
    }

    void removeGrain(size_t grainId)
    {
        // delete grain (swap last)
        grains_[grainId] = grains_[--nGrains_];
        if (nGrains_ == 0 || state_ == kWaitingForLastGrains) {
            state_ = kOff;
        }
    }

    void newGrain()
    {
        auto rate  = synthRate_;
        auto phase = -T(1) + synthPhase_;
        auto delay = dsp::max(T(0), 2.f * (playRate_ - T(1)) / rate);

        addGrain({phase, rate, delay});
    }

    void start()
    {
        if (state_ == kOff) {
            // reset phase
            synthPhase_ = T(1);
        }
        state_ = kOn;
    }

    void stop() { state_ = kWaitingForLastGrains; }

    template <class Ctxt, class DL> void process(Ctxt ctxt, DL delayline)
    {
        delayline.write(ctxt, ctxt.getInput());

        if (state_ == kOff) {
            ctxt.setOutput(T(0));
            return;
        }

        if (synthPhase_ >= T(1)) {
            synthPhase_ -= T(1);
            if (state_ == kOn) {
                newGrain();
            }
        }

        auto x = T(0);
        for (size_t gId = 0; gId < nGrains_; ++gId) {
            auto &grain = grains_[gId];

            if (grain.phase_ > T(1)) {
                removeGrain(gId);
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

            grain.phase_ += grain.rate_;
            grain.delay_ += T(1) - playRate_;
        }

        ctxt.setOutput(x);

        // increment synth phase
        synthPhase_ += synthRate_;
    }

  private:
    State state_{kOff};

    T synthRate_{1.f};
    T synthPhase_{};

    T playRate_{1};

    size_t nGrains_{};
    std::array<Grain, MaxGrains> grains_{};
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
