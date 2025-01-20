#pragma once

#include "Context.h"
#include "Delay.h"
#include "FastMath.h"
#include "Signal.h"
#include <cmath>

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
    static constexpr auto Pad    = Signal<N>::VectorSize;
    static constexpr auto NCoeff = Order + 1;

    template <int N_ = N>
    class DL : public CopyDelayLine<N_, nextAlignedOffset(Order, Pad)>
    {
    };

    static constexpr auto PaddedLength = NCoeff + Pad * 2 - 1;

    FIRFilter() = default;
    FIRFilter(std::array<Signal<N>, NCoeff> &&b)
    {
        for (int k = 0; k < NCoeff; ++k) {
            for (int i = 0; i < N; ++i) {
                b_[PaddedLength - Pad - k][i] = b[k][i];
            }
        }
    }

    template <class Ctxt, class DL>
    __attribute__((always_inline)) void process(Ctxt c, DL &delayline) const
    {
        auto &x = c.getIn();

        typename Ctxt::Type sums[Ctxt::VecSize] = {0};

        for (int j = 0; j < NCoeff + Ctxt::VecSize - 1; ++j) {
            auto n        = j % Ctxt::VecSize;
            auto delay    = j - n;
            auto &x0      = delay == 0 ? x : delayline.read(c, delay);
            const auto &b = b_[PaddedLength - Pad - j].toVector();

#pragma omp simd
            for (int k = 0; k < x0.size(); ++k) {
                for (int i = 0; i < x0[0].size(); ++i) {
                    sums[n][k][i] += x0[k][i] * b[k][i % N];
                }
            }
        }

        delayline.write(c, x);

#pragma omp simd
        for (int n = 0; n < x.size(); ++n) {
            for (int i = 0; i < x[0].size(); ++i) {
                x[n][i] = 0;
                for (int k = 0; k < x.size(); ++k) {
                    x[n][i] += sums[n][k][i];
                }
            }
        }
    }

  private:
    Signal<N> b_[PaddedLength] = {0};
};

template <int N, int Order, int M> class FIRDecimate
{
  public:
    static constexpr auto Pad    = Signal<N>::VectorSize;
    static constexpr auto NCoeff = (Order + 1) * M;

    template <int Offset = 0>
    class DL : public DelayLine<nextAlignedOffset(NCoeff - 1, Pad), Offset>
    {
    };

    static constexpr auto PaddedLength = NCoeff + Pad * 2 - 1;

    FIRDecimate()
    {
        for (int n = 0; n < NCoeff; ++n) {
            for (int i = 0; i < N; ++i) {
                auto mid = NCoeff / 2.f;
                b_[PaddedLength - Pad - n][i] =
                    sinc((n - mid) / (mid)) * sinc((n - mid) / M) / M;
            }
        }
    }

    template <class CtxtIn, class CtxtOut, class DL>
    __attribute__((always_inline)) void decimate(CtxtIn cin, CtxtOut cout,
                                                 DL &delayline) // const
    {
        static_assert(CtxtOut::VecSize == 1);
        auto decimateId = decimateId_;

        contextFor(cin)
        {
            auto x = c.getIn();

            auto xOffset = decimateId;
            while (xOffset < CtxtIn::VecSize) {
                typename CtxtIn::Type sum = {0};

                for (int delay = 0; delay < NCoeff + CtxtIn::VecSize - 1;
                     delay += CtxtIn::VecSize) {
                    auto &x0 = delay == 0 ? x : delayline.read(c, delay);
                    const auto &b =
                        b_[PaddedLength - Pad - delay - xOffset].toVector();

#pragma omp simd
                    for (int k = 0; k < x0.size(); ++k) {
                        for (int i = 0; i < x0[0].size(); ++i) {
                            sum[k][i] += x0[k][i] * b[k][i % N];
                        }
                    }
                }

                auto &xout = cout.getIn();
#pragma omp simd
                for (int i = 0; i < sum[0].size(); ++i) {
                    xout[i] = 0;
                    for (int k = 0; k < sum.size(); ++k) {
                        xout[0][i] += sum[k][i];
                    }
                }

                xOffset += M;
                cout.next();
            }

            decimateId = (decimateId + CtxtIn::VecSize) % M;
            delayline.write(c, x);
        }

        decimateId_ = decimateId;
    }

    // private:
    Signal<N> b_[PaddedLength] = {0};
    int decimateId_            = 0;
};

template <int N, int Order, int M> class FIRInterpolate
{
  public:
    static constexpr auto Pad    = Signal<N>::VectorSize;
    static constexpr auto NCoeff = (Order + 1) * M;

    template <int N_ = N>
    class DL : public CopyDelayLine<N_, nextAlignedOffset(NCoeff - 1, Pad)>
    {
    };

    static constexpr auto PaddedLength = NCoeff + Pad * 2 - 1;

    FIRInterpolate()
    {
        for (int n = 0; n < NCoeff; ++n) {
            for (int i = 0; i < N; ++i) {
                auto mid = NCoeff / 2.f;
                b_[PaddedLength - Pad - n][i] =
                    sinc((n - mid) / (mid)) * sinc((n - mid) / M) / M;
            }
        }
    }

    template <class CtxtIn, class CtxtOut, class DL>
    __attribute__((always_inline)) void decimate(CtxtIn cin, CtxtOut cout,
                                                 DL &delayline) // const
    {
        static_assert(CtxtOut::VecSize == 1);
        auto decimateId = decimateId_;

        contextFor(cin)
        {
            auto x = c.getIn();

            auto xOffset = decimateId;
            while (xOffset < CtxtIn::VecSize) {
                typename CtxtIn::Type sum = {0};

                for (int delay = 0; delay < NCoeff + CtxtIn::VecSize - 1;
                     delay += CtxtIn::VecSize) {
                    auto &x0 = delay == 0 ? x : delayline.read(c, delay);
                    const auto &b =
                        b_[PaddedLength - Pad - delay - xOffset].toVector();

#pragma omp simd
                    for (int k = 0; k < x0.size(); ++k) {
                        for (int i = 0; i < x0[0].size(); ++i) {
                            sum[k][i] += x0[k][i] * b[k][i % N];
                        }
                    }
                }

                auto &xout = cout.getIn();
#pragma omp simd
                for (int i = 0; i < sum[0].size(); ++i) {
                    xout[i] = 0;
                    for (int k = 0; k < sum.size(); ++k) {
                        xout[0][i] += sum[k][i];
                    }
                }

                xOffset += M;
                cout.next();
            }

            decimateId = (decimateId + CtxtIn::VecSize) % M;
            delayline.write(c, x);
        }

        decimateId_ = decimateId;
    }

    // private:
    Signal<N> b_[PaddedLength] = {0};
    int decimateId_            = 0;
};

} // namespace dsp
