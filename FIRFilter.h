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

template <int N, int Order> class FIRFilter
{
  public:
    static constexpr auto Pad    = fSample<N>::VectorSize;
    static constexpr auto NCoeff = Order + 1;

    template <int N_ = N>
    class DL : public CopyDelayLine<N_, nextAlignedOffset(Order, Pad)>
    {
    };

    static constexpr auto PaddedLength = NCoeff + Pad * 2 - 1;

    FIRFilter() = default;
    FIRFilter(const std::array<fData<N>, NCoeff> b)
    {
        for (size_t k = 0; k < NCoeff; ++k) {
            for (size_t i = 0; i < N; ++i) {
                b_[PaddedLength - Pad - k][i] = b[k][i];
            }
        }
    }

    template <class Ctxt, class DL>
    __attribute__((always_inline)) void process(Ctxt c, DL &delayline) const
    {
        auto &x = c.getSignal();

        typename Ctxt::Type sums[Ctxt::VecSize] = {0};

        for (auto j = 0; j < NCoeff + Ctxt::VecSize - 1; ++j) {
            auto n        = j % Ctxt::VecSize;
            auto delay    = j - n;
            auto &x0      = delay == 0 ? x : delayline.read(c, delay);
            const auto &b = b_[PaddedLength - Pad - j].toVector();

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
    fSample<N> b_[PaddedLength] = {0};
};

template <int N, int Order, int M> class FIRDecimate
{
  public:
    static constexpr auto Pad    = fSample<N>::VectorSize;
    static constexpr auto NCoeff = (Order + 1) * M;

    template <int Offset = 0>
    class DL : public DelayLine<nextAlignedOffset(NCoeff - 1, Pad), Offset>
    {
    };

    static constexpr auto PaddedLength = NCoeff + Pad * 2 - 1;

    FIRDecimate(float cutoff = 1.0f)
    {
        for (auto n = 0; n < NCoeff; ++n) {
            for (auto i = 0; i < N; ++i) {
                auto mid  = NCoeff / 2.f;
                auto freq = cutoff / M;
                b_[PaddedLength - Pad - n][i] =
                    window::Kaiser<140>::generate((n - mid) / (mid)) *
                    sinc((n - mid) * freq) * freq;
            }
        }
    }

    template <class CtxtIn, class CtxtOut, class DL>
    __attribute__((always_inline)) int
    decimate(CtxtIn cin, CtxtOut &cout, DL &delayline, int decimateId) const
    {
        static_assert(CtxtOut::VecSize == 1);
        auto decimatedBlockSize = 0;
        auto id                 = decimateId;
        auto coutCopy           = cout;

        contextFor(cin)
        {
            auto x = c.getSignal();

            size_t xOffset = (M - id) % M;
            while (xOffset < CtxtIn::VecSize) {
                typename CtxtIn::Type sum = {0};

                for (size_t delay = 0; delay < NCoeff + CtxtIn::VecSize - 1;
                     delay += CtxtIn::VecSize) {
                    auto &x0 = delay == 0 ? x : delayline.read(c, delay);
                    const auto &b =
                        b_[PaddedLength - Pad - delay - xOffset].toVector();

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

            id = (id + CtxtIn::VecSize) % M;
            delayline.write(c, x);
        }

        cout.setBlockSize(decimatedBlockSize);
        return id;
    }

  private:
    fSample<N> b_[PaddedLength] = {0};
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
template <int N, int Order, int L> class FIRInterpolate
{
  public:
    static constexpr auto Pad    = fSample<N>::VectorSize;
    static constexpr auto NCoeff = (Order + 1);

    template <int N_ = N>
    class DL : public CopyDelayLine<N_, nextAlignedOffset(NCoeff, Pad)>
    {
    };

    static constexpr auto PaddedLength = NCoeff + Pad * 2 - 1;

    FIRInterpolate(float cutoff = 1.f)
    {
        for (size_t l = 0; l < L; ++l) {
            for (size_t n = 0; n < NCoeff; ++n) {
                for (size_t i = 0; i < N; ++i) {
                    auto freq = cutoff / L;
                    auto mid  = NCoeff * L / 2.f;
                    auto k    = l + n * L;
                    b_[l][PaddedLength - Pad - n][i] =
                        window::Kaiser<140>::generate((k - mid) / (mid)) *
                        sinc((k - mid) * freq) * cutoff;
                }
            }
        }
    }

    template <class CtxtIn, class CtxtOut, class DL>
    __attribute__((always_inline)) int interpolate(CtxtIn cin, CtxtOut cout,
                                                   DL &delayline,
                                                   int interpolateId) const
    {
        static_assert(CtxtIn::VecSize == 1);
        static_assert(CtxtOut::VecSize == 1);
        auto id = interpolateId;

        contextFor(cout)
        {
            if (id == 0) {
                auto x = cin.getSignal();
                delayline.write(cin, x);
                cin.next();
            }

            typename CtxtIn::Type sum = {0};

            for (size_t delay = 0; delay < NCoeff + CtxtIn::VecSize - 1;
                 delay += CtxtIn::VecSize) {
                // if delayline change to external buffer delayline we will have
                // a problem because CopyDelayLine index changes after each
                // write whereas Buffer Delay Line changes with new context...
                // we will see in the future
                auto &x0      = delayline.read(c, delay + 1);
                const auto &b = b_[id][PaddedLength - Pad - delay].toVector();

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
        return id;
    }

  private:
    fSample<N> b_[L][PaddedLength] = {0};
};

} // namespace dsp
