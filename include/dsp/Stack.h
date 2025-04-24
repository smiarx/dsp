#pragma once

#include "Signal.h"
#include "Utils.h"

namespace dsp
{

template <size_t N, size_t Size> class Stack
{
    static_assert(isPow2(Size), "Size must be a power of 2");
    static constexpr auto kMask = Size - 1;

  public:
    fSample<N> push(const fSample<N> &x)
    {
        n_         = (n_ + 1) & kMask;
        auto y     = stack_[n_];
        stack_[n_] = x;
        return y;
    }

    [[nodiscard]] fSample<N> get() const { return stack_[n_]; }

    [[nodiscard]] const fSample<N> *getSamples() const { return stack_; }
    [[nodiscard]] const size_t *getPos() const { return &n_; }

  private:
    size_t n_{Size - 1};
    fSample<N> stack_[Size]{};
};

} // namespace dsp
