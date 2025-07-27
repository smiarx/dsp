#pragma once

#include "Context.h"

namespace dsp
{
template <typename T> class DoubleRamp
{
  public:
    enum State {
        kOff = 0,
        kUp,
        kDown,
    };

    void set(const T &maxValue, const T &up, const T &down)
    {
        target_   = maxValue;
        stepUp_   = (maxValue - value_) / up;
        stepDown_ = (-load(maxValue)) / down;
        state_    = kUp;
    }

    bool isRunning() { return !all(load(state_) == kOff); }

    auto process()
    {
        auto value = load(value_);

        if (isRunning()) {
            auto target = load(target_);
            auto state  = load(state_);

            auto isStateUp   = state == kUp;
            auto isStateDown = load(state_) == kDown;

            // up part
            auto valueUp = value + stepUp_;

#ifdef DSP_ARM32
            auto greaterThanTarget = abs(valueUp) + 1e-6f >= abs(target);
#else
            auto greaterThanTarget = abs(valueUp) >= abs(target);
#endif
            valueUp = blend(greaterThanTarget, target * 2 - valueUp, valueUp);
            state   = blend(greaterThanTarget, load(intType<T>(kDown)), state);

            value = blend(isStateUp, valueUp, value);

            // down part
            auto valueDown = value + stepDown_;

            auto changedSign =
                signbit(valueDown) ^ signbit(target) || abs(valueDown) < 1e-5;
            valueDown = blend(changedSign, load(T(0)), valueDown);
            state     = blend(changedSign, load(intType<T>(kOff)), state);

            value = blend(isStateDown, valueDown, value);

            value_ = value;
            state_ = state;
        }

        return value;
    }

  private:
    T value_{};
    T target_{};
    T stepUp_{};
    T stepDown_{};
    intType<T> state_{};
};
} // namespace dsp
