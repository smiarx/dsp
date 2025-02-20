#pragma once

#include "Context.h"
#include "Signal.h"
#include <cassert>
#include <cmath>

namespace dsp
{

template <int N = 1, bool useVector = false> class ControlSmoother
{
    /* smooth control values between audio blocks */
    using Type    = Signal<N>;
    using OutType = typename std::conditional<useVector, typename Type::Vector,
                                              typename Type::Scalar>::type;

    /* linear smoother for control values */
  public:
    ControlSmoother(Type target) : target_{target}
    {
        arrayFor(value_, k) { value_[k] = target; }
    }

    void set(Type target, float invBlockSize)
    {
        static constexpr auto VecSize = OutType().size();
        if (target_ != target) {
            target_ = target;
#pragma omp simd
            inFor(value_, k, i)
            {
                auto step = (target_[i] - value_[k][i]) * invBlockSize;
                value_[k][i] -= step * (VecSize - 1 - k);
                step_[k][i] = step * VecSize;
            }
            active_ = true;
        }
    }

    auto get() const { return value_; }

    bool isActive() { return active_; }

    bool step()
    {
        if (!active_) return false;
        inFor(value_, k, i) { value_[k][i] += step_[k][i]; }
        return true;
    }

    void reset()
    {
        active_ = false;
        arrayFor(value_, k) { value_[k] = target_; }
    }

  private:
    Type target_{0.f};
    OutType value_{0.f};
    OutType step_{0.f};
    bool active_{false};
};

class SmootherExp
{
  public:
    SmootherExp(float target) : target_{target}, value_(target) {}
    void set(float target) { target_ = target; }
    void setTime(float time) { coef_ = 1.f - powf(0.001f, 1.f / time); }
    float step()
    {
        value_ += coef_ * (target_ - value_);
        return value_;
    }

  private:
    float target_{0.f};
    float value_{0.f};
    float coef_{0.f};
};

template <int N> class SmootherLin
{
  public:
    SmootherLin() = default;
    SmootherLin(Signal<N> target) : value_{target} {}

    void set(Signal<N> target, int count)
    {
        for (int i = 0; i < N; ++i) {
            step_[i] = (target[i] - value_[i]) / count;
            count_   = count;
        }
    }

    Signal<N> step()
    {
        if (count_ == 0) {
            step_ = {0.f};
        } else {
            for (int i = 0; i < N; ++i) {
                value_[i] += step_[i];
            }
            --count_;
        }
        return value_;
    }

  private:
    Signal<N> value_{{0.f}};
    Signal<N> step_{{0.f}};
    int count_{0};
};

} // namespace dsp
