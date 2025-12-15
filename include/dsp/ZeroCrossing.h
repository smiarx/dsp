#pragma once

#include "dsp/simd/multi.h"
#include "dsp/simd/simd_math.h"

namespace dsp
{
template <typename T> class ZeroCrossing
{
  public:
    template <class Ctxt> auto process(const Ctxt ctxt)
    {
        auto x    = ctxt.getInput();
        auto sign = dsp::signbit(x);

        auto result = sign ^ prevsign_;
        prevsign_   = sign;

        return result;
    }

  private:
    decltype(dsp::signbit(dsp::load(T()))) prevsign_{};
};
} // namespace dsp
