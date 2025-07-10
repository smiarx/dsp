#pragma once

#include "Context.h"
#include "MultiVal.h"
#include <cassert>
#include <cmath>

namespace dsp
{

template <typename T, bool Vectorize = false> class ControlSmoother
{
    /* smooth control values between audio blocks */
    using outT                   = std::conditional_t<Vectorize, batch<T>, T>;
    using bT                     = baseType<T>;
    static constexpr auto kWidth = kTypeWidth<outT> / kTypeWidth<T>;

    /* linear smoother for control values */
  public:
    ControlSmoother() = default;
    ControlSmoother(const T &target) : target_{target}, value_(target) {}

    void set(const T &target, bT invBlockSize)
    {
        if (!all(load(target_) == load(target))) {
            target_ = target;

            auto step = (target - load(value_)) * invBlockSize;
            step_[0]  = step;

            // create incr vector (eg. {1,2,3,4})
            T incr[kWidth];
            for (size_t i = 0; i < kWidth; ++i) incr[i] = i + 1;

            Context<T, Vectorize> ctxt(step_.data());
            ctxt.setOutput(ctxt.load(*incr) * step_[0]);

            active_ = true;
        }
    }

    [[nodiscard]] auto getTarget() const { return load(target_); }

    bool isActive() { return active_; }

    template <class Ctxt> auto step(const Ctxt &)
    {
        static_assert(
            Vectorize || !Ctxt::kUseVec,
            "Cannot use scalar ControlSmoother in Vectorized context");

        Context<T, Ctxt::kUseVec> c{nullptr};

        using sigT               = decltype(c.getInput());
        constexpr auto kIncrSize = Ctxt::kIncrSize;

        // trick to broadcast value to correct type
        auto broadcastValue = sigT(0) + value_;

        // store value in array
        T values[kIncrSize];
        c.setData(values);
        c.setOutput(broadcastValue);

        if (active_) {
            // add step if active
            c.setOutput(c.getInput() + c.load(*step_.data()));
        }

        value_ = values[kIncrSize - 1];

        return c.load(*values);
    }

    [[nodiscard]] auto get() { return load(value_); }

    void reset()
    {
        active_ = false;
        value_  = target_;
    }

  private:
    T target_{};
    T value_{};

    std::array<T, kWidth> step_{};

    bool active_{false};
};

template <typename T> class SmootherExp
{
    using bt = baseType<T>;

  public:
    SmootherExp(const T &target) : target_{target}, value_(target) {}
    void set(const T &target) { target_ = target; }
    void setTime(bt time) { coef_ = bt(1) - std::pow(bt(0.001), bt(1) / time); }
    auto step()
    {
        value_ += coef_ * (target_ - value_);
        return load(value_);
    }

  private:
    T target_{};
    T value_{};
    bt coef_{};
};

template <typename T> class SmootherLin
{
  public:
    SmootherLin() = default;
    SmootherLin(const T &target) : value_{target} {}

    void set(const T &target, int count)
    {
        step_  = (load(target) - value_) / count;
        count_ = count;
    }

    auto step()
    {
        if (count_ == 0) {
            step_ = {};
        } else {
            value_ += load(step_);
            --count_;
        }
        return load(value_);
    }

  private:
    T value_{};
    T step_{};
    int count_{};
};

} // namespace dsp
