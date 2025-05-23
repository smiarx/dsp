#pragma once

#include "Context.h"
#include "Delay.h"
#include "FastMath.h"
#include "Signal.h"
#include "Window.h"
#include <cmath>
#include <cstdio>

/*
    x0 |x1 x2 x3 x4| x5 x6 x7 x8 |x9
       |           |             |
 b6 b5 |b4 b3 b2 b1| b0 0  0  0  |
       |           |             |
 b7 b6 |b5 b4 b3 b2| b1 b0 0  0  |
       |           |             |
    b7 |b6 b5 b4 b3| b2 b1 b0 0  |
       |           |             |
       |b7 b6 b5 b4| b3 b2 b1 b0 |

*/

namespace dsp
{

template <size_t N, size_t Order> class FIRFilter
{
  public:
    static constexpr auto kPad    = fSample<N>::kVectorSize;
    static constexpr auto kNCoeff = Order + 1;

    template <size_t N_ = N>
    class DL : public CopyDelayLine<N_, nextAlignedOffset(Order, kPad)>
    {
    };

    static constexpr auto kPaddedLength = kNCoeff + kPad * 2 - 1;

    FIRFilter() = default;
    FIRFilter(const std::array<fData<N>, kNCoeff> b)
    {
        for (size_t k = 0; k < kNCoeff; ++k) {
            for (size_t i = 0; i < N; ++i) {
                b_[kPaddedLength - kPad - k][i] = b[k][i];
            }
        }
    }

    template <class Ctxt, class DL>
    __attribute__((always_inline)) void process(Ctxt c, DL &delayline) const
    {
        auto &x = c.getSignal();

        typename Ctxt::Type sums[Ctxt::kVecSize] = {};

        for (size_t j = 0; j < kNCoeff + Ctxt::kVecSize - 1; ++j) {
            auto n        = j % Ctxt::kVecSize;
            auto delay    = j - n;
            auto &x0      = delay == 0 ? x : delayline.read(c, delay);
            const auto &b = b_[kPaddedLength - kPad - j].toVector();

#pragma omp simd
            for (size_t k = 0; k < x0.size(); ++k) {
                for (size_t i = 0; i < x0[0].size(); ++i) {
                    sums[n][k][i] += x0[k][i] * b[k][i % N];
                }
            }
        }

        delayline.write(c, x);

#pragma omp simd
        for (size_t n = 0; n < x.size(); ++n) {
            for (size_t i = 0; i < x[0].size(); ++i) {
                x[n][i] = 0;
                for (size_t k = 0; k < x.size(); ++k) {
                    x[n][i] += sums[n][k][i];
                }
            }
        }
    }

  private:
    fSample<N> b_[kPaddedLength] = {};
};

template <size_t N, size_t Order, size_t M> class FIRDecimate
{
  public:
    static constexpr auto kPad    = fSample<N>::kVectorSize;
    static constexpr auto kNCoeff = (Order + 1) * M;

    template <int Offset = 0>
    class DL : public DelayLine<nextAlignedOffset(kNCoeff - 1, kPad), Offset>
    {
    };

    static constexpr auto kPaddedLength = kNCoeff + kPad * 2 - 1;

    FIRDecimate(float cutoff = 1.0f)
    {
        for (size_t n = 0; n < kNCoeff; ++n) {
            for (size_t i = 0; i < N; ++i) {
                auto mid  = kNCoeff / 2.f;
                auto freq = cutoff / M;
                auto fn   = static_cast<float>(n);
                b_[kPaddedLength - kPad - n][i] =
                    window::Kaiser<140>::generate((fn - mid) / (mid)) *
                    sinc((fn - mid) * freq) * freq;
            }
        }
    }

    template <class CtxtIn, class CtxtOut, class DL>
    __attribute__((always_inline)) int
    decimate(CtxtIn cin, CtxtOut &cout, DL &delayline, int decimateId) const
    {
        static_assert(CtxtOut::kVecSize == 1);
        auto decimatedBlockSize = 0;
        auto id                 = static_cast<size_t>(decimateId);
        auto coutCopy           = cout;

        contextFor(cin)
        {
            auto x = c.getSignal();

            size_t xOffset = (M - id) % M;
            while (xOffset < CtxtIn::kVecSize) {
                typename CtxtIn::Type sum = {};

                for (size_t delay = 0; delay < kNCoeff + CtxtIn::kVecSize - 1;
                     delay += CtxtIn::kVecSize) {
                    auto &x0 = delay == 0
                                   ? x
                                   : delayline.read(c, static_cast<int>(delay));
                    const auto &b =
                        b_[kPaddedLength - kPad - delay - xOffset].toVector();

#pragma omp simd
                    for (size_t k = 0; k < x0.size(); ++k) {
                        for (size_t i = 0; i < x0[0].size(); ++i) {
                            sum[k][i] += x0[k][i] * b[k][i % N];
                        }
                    }
                }

                auto &xout = coutCopy.getSignal();
#pragma omp simd
                for (size_t i = 0; i < sum[0].size(); ++i) {
                    xout[0][i] = 0;
                    for (size_t k = 0; k < sum.size(); ++k) {
                        xout[0][i] += sum[k][i];
                    }
                }

                ++decimatedBlockSize;
                xOffset += M;
                coutCopy.next();
            }

            id = (id + CtxtIn::kVecSize) % M;
            delayline.write(c, x);
        }

        cout.setBlockSize(decimatedBlockSize);
        return static_cast<int>(id);
    }

  private:
    fSample<N> b_[kPaddedLength] = {};
};

/*
    x0 |x1 x2 x3 x4| x5 x6 x7 x8 |x9
    0  |0  0  b6 b3| b0 0  0  0  |
    0  |0  0  b7 b4| b1 0  0  0  |
    0  |0  0  b8 b5| b2 0  0  0  |
    0  |0  0  0  b6| b3 0  0  0  |
    0  |0  0  0  b7| b4 0  0  0  |
    0  |0  0  0  b8| b5 0  0  0  |
    0  |0  0  0  0 | b6 0  0  0  |
    0  |0  0  0  0 | b7 0  0  0  |


 */
template <size_t N, size_t Order, size_t L> class FIRInterpolate
{
  public:
    static constexpr auto kPad    = fSample<N>::kVectorSize;
    static constexpr auto kNCoeff = (Order + 1);

    template <size_t N_ = N>
    class DL : public CopyDelayLine<N_, nextAlignedOffset(kNCoeff, kPad)>
    {
    };

    static constexpr auto kPaddedLength = kNCoeff + kPad * 2 - 1;

    FIRInterpolate(float cutoff = 1.f)
    {
        for (size_t l = 0; l < L; ++l) {
            for (size_t n = 0; n < kNCoeff; ++n) {
                for (size_t i = 0; i < N; ++i) {
                    auto freq = cutoff / L;
                    auto mid  = kNCoeff * L / 2.f;
                    auto k    = l + n * L;
                    auto fk   = static_cast<float>(k);
                    b_[l][kPaddedLength - kPad - n][i] =
                        window::Kaiser<140>::generate((fk - mid) / (mid)) *
                        sinc((fk - mid) * freq) * cutoff;
                }
            }
        }
    }

    template <class CtxtIn, class CtxtOut, class DL>
    __attribute__((always_inline)) int interpolate(CtxtIn cin, CtxtOut cout,
                                                   DL &delayline,
                                                   int interpolateId) const
    {
        static_assert(CtxtIn::kVecSize == 1);
        static_assert(CtxtOut::kVecSize == 1);
        auto id = static_cast<size_t>(interpolateId);

        contextFor(cout)
        {
            if (id == 0) {
                auto x = cin.getSignal();
                delayline.write(cin, x);
                cin.next();
            }

            typename CtxtIn::Type sum = {};

            for (size_t delay = 0; delay < kNCoeff + CtxtIn::kVecSize - 1;
                 delay += CtxtIn::kVecSize) {
                // if delayline change to external buffer delayline we will have
                // a problem because CopyDelayLine index changes after each
                // write whereas Buffer Delay Line changes with new context...
                // we will see in the future
                auto &x0      = delayline.read(c, static_cast<int>(delay + 1));
                const auto &b = b_[id][kPaddedLength - kPad - delay].toVector();

#pragma omp simd
                for (size_t k = 0; k < x0.size(); ++k) {
                    for (size_t i = 0; i < x0[0].size(); ++i) {
                        // y_n = b0 * x_n + b1 * x_{n-1} + b2 *x_{n-2} + ....
                        sum[k][i] += x0[k][i] * b[k][i % N];
                    }
                }
            }

            auto &xout = c.getSignal();
#pragma omp simd
            for (size_t i = 0; i < sum[0].size(); ++i) {
                xout[0][i] = 0.f;
                for (size_t k = 0; k < sum.size(); ++k) {
                    xout[0][i] += sum[k][i];
                }
            }

            id = (id + 1) % L;
        }
        return static_cast<int>(id);
    }

  private:
    fSample<N> b_[L][kPaddedLength] = {};
};

} // namespace dsp
