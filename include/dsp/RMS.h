#pragma once

#include "Delay.h"
#include "Signal.h"

namespace dsp
{

template <size_t N, size_t Size, size_t Overlap = 0> class RMS
{
  public:
    static constexpr auto Shift           = Size - Overlap;
    static constexpr bool ShiftDivideSize = Size % Shift == 0;
    static constexpr auto NOverlaps = ShiftDivideSize ? Size / Shift : Size;

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

                if constexpr (!ShiftDivideSize) {
                    moveMeanSq(c);
                }

                n_ = (n_ + 1) % Shift;
                if (n_ == 0) {
                    if constexpr (ShiftDivideSize) {
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
        static_assert(Ctxt::VecSize == 1);

        auto &lastOverlap = overlaps_.tail(c)[0];

        arrayFor(sumsq_, i) { sumsq_[i] += firstOverlap_[i] - lastOverlap[i]; }

        overlaps_.write(c, firstOverlap_.toSignal());
        firstOverlap_ = {};
    }

    int n_ = 0;
    fSample<N> firstOverlap_{};
    CopyDelayLine<N, NOverlaps> overlaps_{};
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
    int n_ = 0;
    fSample<N> sumsq_{};
};
} // namespace dsp
