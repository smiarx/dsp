#pragma once

#include "Delay.h"

#ifdef __GNUC__
#ifndef __clang__
// gcc fast-math has unpredictible behavior here
#pragma GCC push_options
#pragma GCC optimize("no-fast-math")
#endif
#endif

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T, size_t Size, size_t Overlap = 0> class RMS
{
  public:
    static constexpr size_t kShift         = Size - Overlap;
    static constexpr bool kShiftDivideSize = Size % kShift == 0;
    static constexpr auto kNOverlaps = kShiftDivideSize ? Size / kShift : Size;

    template <class Ctxt, class Stack>
    void processBlock(Ctxt ctxt, Stack &stack)
    {
        CTXTRUN(ctxt)
        {
            auto x = ctxt.getInput();

            firstOverlap_ += x * x;
            if constexpr (!kShiftDivideSize) {
                moveMeanSq(ctxt);
            }

            n_ = (n_ + 1) % kShift;
            if (n_ == 0) {
                if constexpr (kShiftDivideSize) {
                    moveMeanSq(ctxt);
                }
                auto rms = sqrt(load(sumsq_) / static_cast<baseType<T>>(Size));
                stack.push(rms);
            }
        };
    }

  private:
    template <class Ctxt> void moveMeanSq(Ctxt c)
    {
        static_assert(!Ctxt::kUseVec);

        auto lastOverlap = overlaps_.tail(c);

        sumsq_ = max(load(T(0)), load(sumsq_) + firstOverlap_ - lastOverlap);

        overlaps_.write(c, firstOverlap_);
        firstOverlap_ = {};
    }

    size_t n_ = 0;
    T firstOverlap_{};
    CopyDelayLine<T, kNOverlaps> overlaps_{};
    T sumsq_{};
};

/* no overlap */
template <typename T, size_t Size> class RMS<T, Size>
{
  public:
    template <class Ctxt, class Stack>
    void processBlock(Ctxt ctxt, Stack &stack)
    {
        CTXTRUN(ctxt)
        {
            auto x = ctxt.getInput();

            sumsq_ += x * x;
            n_ = (n_ + 1) % Size;

            if (n_ == 0) {
                auto rms = sqrt(load(sumsq_) / static_cast<baseType<T>>(Size));
                stack.push(rms);
                sumsq_ = {};
            }
        };
    }

  private:
    size_t n_ = 0;
    T sumsq_{};
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp

#ifdef __GNUC__
#ifndef __clang__
#pragma GCC pop_options
#endif
#endif
