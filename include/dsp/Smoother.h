#pragma once

#include "Context.h"
#include "Signal.h"
#include <cassert>
#include <cmath>

namespace dsp
{

template <size_t N = 1, bool Vectorize = false> class ControlSmoother
{
    /* smooth control values between audio blocks */
    using Type    = fSample<N>;
    using OutType = std::conditional_t<Vectorize, typename Type::Vector,
                                       typename Type::Scalar>;

    /* linear smoother for control values */
  public:
    ControlSmoother(Type target) : target_{target}
    {
        arrayFor(value_, k) { value_[k] = target; }
    }

    void set(Type target, float invBlockSize)
    {
        static constexpr auto kVecSize = OutType().size();
        if (target_ != target) {
            target_ = target;
#pragma omp simd
            inFor(value_, k, i)
            {
                auto step = (target_[i] - value_[k][i]) * invBlockSize;
                value_[k][i] -= step * static_cast<float>(kVecSize - 1 - k);
                step_[k][i] = step * kVecSize;
            }
            active_ = true;
        }
    }

    [[nodiscard]] auto get() const { return value_; }
    [[nodiscard]] auto getTarget() const { return target_; }

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
    OutType value_{};
    OutType step_{};
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

template <size_t N> class SmootherLin
{
  public:
    SmootherLin() = default;
    SmootherLin(fData<N> target) : value_{target} {}

    void set(fData<N> target, int count)
    {
        auto fcount = static_cast<float>(count);
        for (size_t i = 0; i < N; ++i) {
            step_[i] = (target[i] - value_[i]) / fcount;
            count_   = count;
        }
    }

    fData<N> step()
    {
        if (count_ == 0) {
            step_ = {0.f};
        } else {
            for (size_t i = 0; i < N; ++i) {
                value_[i] += step_[i];
            }
            --count_;
        }
        return value_;
    }

  private:
    fData<N> value_{{0.f}};
    fData<N> step_{{0.f}};
    int count_{0};
};

} // namespace dsp
