#pragma once

#include "Context.h"
#include "Signal.h"
#include <cassert>
#include <cmath>

namespace dsp
{

class Smoother
{
  public:
    Smoother(float target) : target_{target}, value_(target) {}
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

// class ControlSmoother
//{
//   public:
//     ControlSmoother(float target) : smoother_{target},
//     linearSmoother_{target}
//     {
//     }
//     void set(float target) { smoother_.set(target); }
//     void setTime(float time)
//     {
//         smoother_.setTime(time);
//     }
//     void controlStep(float invBlockSize)
//     {
//         linearSmoother_.set(smoother_.step(), invBlockSize);
//     }
//     bool audioStep() { return linearSmoother_.step(); }
//     float get() const { return linearSmoother_.get(); }
//     void reset() { linearSmoother_.reset(); }
//
//   private:
//     Smoother smoother_;
//     LinearSmoother linearSmoother_;
// };

} // namespace dsp
