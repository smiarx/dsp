#pragma once

#include "MultiVal.h"
#include "Utils.h"

namespace dsp
{

template <typename T, size_t Size> class Stack
{
    static_assert(isPow2(Size), "Size must be a power of 2");
    static constexpr auto kMask = Size - 1;

  public:
    auto push(const T &x)
    {
        n_         = (n_ + 1) & kMask;
        auto y     = stack_[n_];
        stack_[n_] = x;
        return y;
    }

    [[nodiscard]] auto get() const { return loadfuncs::load(stack_[n_]); }

    [[nodiscard]] const T *getSamples() const { return stack_; }
    [[nodiscard]] const size_t *getPos() const { return &n_; }

  private:
    size_t n_{Size - 1};
    T stack_[Size]{};
};

} // namespace dsp
