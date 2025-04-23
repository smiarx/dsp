#pragma once

#include "Delay.h"
#include "Signal.h"

namespace dsp
{

template <size_t N, size_t Size, size_t Overlap = 0> class RMS
{
  public:
    static constexpr size_t kShift         = Size - Overlap;
    static constexpr bool kShiftDivideSize = Size % kShift == 0;
    static constexpr auto kNOverlaps = kShiftDivideSize ? Size / kShift : Size;

    template <class Ctxt, class Stack>
    void processBlock(Ctxt ctxt, Stack &stack)
    {
        contextFor(ctxt)
        {
            auto &x = c.getSignal();

            arrayFor(x, k)
            {
                auto &xk = x[k];
                arrayFor(xk, i) { firstOverlap_[i] += xk[i] * xk[i]; }

                if constexpr (!kShiftDivideSize) {
                    moveMeanSq(c);
                }

                n_ = (n_ + 1) % kShift;
                if (n_ == 0) {
                    if constexpr (kShiftDivideSize) {
                        moveMeanSq(c);
                    }

                    decltype(sumsq_) rms;
                    arrayFor(sumsq_, i)
                    {
                        rms[i] = sqrtf(sumsq_[i] / float(Size));
                    }
                    stack.push(rms);
                }
            }
        }
    }

  private:
    template <class Ctxt> void moveMeanSq(Ctxt c)
    {
        static_assert(Ctxt::kVecSize == 1);

        auto &lastOverlap = overlaps_.tail(c)[0];

        arrayFor(sumsq_, i)
        {
            sumsq_[i] += firstOverlap_[i] - lastOverlap[i];
            sumsq_[i] = std::max(
                0.f, sumsq_[i]); // protect from really small negative numbers
        }

        overlaps_.write(c, firstOverlap_.toSignal());
        firstOverlap_ = {};
    }

    size_t n_ = 0;
    fSample<N> firstOverlap_{};
    CopyDelayLine<N, kNOverlaps> overlaps_{};
    fSample<N> sumsq_{};
};

/* no overlap */
template <size_t N, size_t Size> class RMS<N, Size>
{
  public:
    template <class Ctxt, class Stack>
    void processBlock(Ctxt ctxt, Stack &stack)
    {
        contextFor(ctxt)
        {
            auto &x = c.getSignal();

            arrayFor(x, k)
            {
                auto &xk = x[k];
                arrayFor(xk, i) { sumsq_[i] += xk[i] * xk[i]; }

                n_ = (n_ + 1) % Size;
                if (n_ == 0) {

                    decltype(sumsq_) rms;
                    arrayFor(sumsq_, i)
                    {
                        rms[i] = sqrtf(sumsq_[i] / float(Size));
                    }
                    stack.push(rms);
                    sumsq_ = {};
                }
            }
        }
    }

  private:
    size_t n_ = 0;
    fSample<N> sumsq_{};
};
} // namespace dsp
