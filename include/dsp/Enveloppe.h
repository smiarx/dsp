#pragma once

#include "Context.h"
#include "Signal.h"
#include <cmath>

namespace dsp
{
template <std::size_t N> class DoubleRamp
{
  public:
    enum State {
        kOff = 0,
        kUp,
        kDown,
    };

    void set(fData<N> maxValue, fData<N> up, fData<N> down)
    {
        arrayFor(maxValue, i)
        {
            target_[i]   = maxValue[i];
            stepUp_[i]   = (maxValue[i] - value_[i]) / up[i];
            stepDown_[i] = -maxValue[i] / down[i];
            state_[i]    = kUp;
        }
    }

    bool isRunning()
    {
        bool running = false;
        arrayFor(state_, i) { running |= (state_[i] != kOff); }
        return running;
    }

    fSample<N> process()
    {
        fSample<N> x{};

        arrayFor(x, i)
        {
            if (state_[i] == kUp) {
                value_[i] += stepUp_[i];
                if (std::abs(value_[i]) >= std::abs(target_[i])) {
                    value_[i] = 2 * target_[i] - value_[i];
                    state_[i] = kDown;
                }
            } else if (state_[i] == kDown) {
                value_[i] += stepDown_[i];
                if (std::signbit(value_[i]) ^ std::signbit(target_[i]) ||
                    std::abs(value_[i]) < 1e-5f) {
                    value_[i] = 0.f;
                    state_[i] = kOff;
                }
            }
            x[i] = value_[i];
        }
        return x;
    }

  private:
    fSample<N> value_{};
    fData<N> target_{};
    fData<N> stepUp_{};
    fData<N> stepDown_{};
    std::array<State, N> state_{};
};
} // namespace dsp
