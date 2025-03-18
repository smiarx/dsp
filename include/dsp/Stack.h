#pragma once

#include "Signal.h"
#include "Utils.h"

namespace dsp
{

template <size_t N, size_t Size> class Stack
{
    static_assert(isPow2(Size), "Size must be a power of 2");
    static constexpr auto Mask = Size - 1;

  public:
    fSample<N> push(fSample<N> x)
    {
        n_         = (n_ + 1) & Mask;
        auto y     = stack_[n_];
        stack_[n_] = x;
        return y;
    }

    fSample<N> get() const { return stack_[n_]; }

    const fSample<N> *getSamples() const { return stack_; }
    const size_t *getPos() const { return &n_; }

  private:
    size_t n_{Size - 1};
    fSample<N> stack_[Size]{};
};

} // namespace dsp
