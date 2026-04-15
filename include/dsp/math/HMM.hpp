#pragma once

#include <dsp/cpu/defines.h>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T> class HMM2
{
    // simple two states HMM
  public:
    HMM2(T transition) : transition_{transition} {}

    T getState() const { return state_; };

    void reset()
    {
        state_ = 0.f;
    }

    void process(T observation1, T observation2)
    {
        auto st  = state_ * transition_;
        auto ps2 = state_ + transition_ - st -
                   st; // probability of being in state2 given previous state
        auto ps1o =
            observation1 * (T(1) - ps2); // probability state1 given observation
        auto ps2o = observation2 * ps2;  // probability state2 given observation

        state_ = ps1o / (ps1o + ps2o);
    }

  private:
    T transition_{}; // state transition from 1 -> 1, other transitions are
                     // implied
    T state_{T(0.5)};
};
} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
