#pragma once

#include "Delay.h"
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
        __attribute__((always_inline))
    void process(Ctxt c, DL &delayline) const
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

// enum class PolyLength {
//     Fixed,
//     Var,
// };
//
// template <int Order, int M, PolyLength FixedVar = PolyLength::Fixed>
// class FIRPolyphase
//{
//   public:
//     template <class Ctxt> class MultiRateContext : public Ctxt
//     {
//         friend FIRPolyphase;
//
//         MultiRateContext(Ctxt c, FIRPolyphase &ply) :
//             Ctxt(c), ply_(ply), multiRateId_(ply_.id)
//         {
//         }
//
//         void next(int incr = 1)
//         {
//             if (multiRateId_ == ply_.getFactor()) {
//                 auto incrNewRate = multiRateId_ / getFactor();
//                 Ctxt::next(incrNewRate);
//                 Ctxt::nextIn(incr - incrNewRate);
//
//                 multiRateId_ %= ply_.getFactor();
//             } else {
//                 Ctxt::nextIn(incr);
//             }
//         }
//
//         int getMultiRateId() const { return multiRateId_; }
//
//         void save() { ply_.id_ = multiRateId_; }
//
//       private:
//         const FIRPolyphase &ply_;
//         int multiRateId_;
//     };
//
//     template <class Ctxt>
//     MultiRateContext<Ctxt> getMultiRateContext(Ctxt c) const
//     {
//         return MultiRateContext(c);
//     }
//
//     template <class Ctxt>
//     Ctxt getContext(Ctxt c) { return Ctxt(c, getFactor()); }
//
//     template <int N> struct delayline {
//         Signal<N> accumulator_;
//         std::array<typename FIRFilter<Order>::template delayline<N>, M>
//         components_;
//     };
//
//     /* construct polyphase filter */
//     FIRPolyphase(int factor, float cutoff = 1.f)
//     {
//         int N = factor*Order;
//         for(int m = 0; m < M; ++m)
//         {
//             Signal<Order> b;
//             for(int k = 0; k < Order; ++k)
//             {
//                 int n = m+k*M;
//                 auto x = n - N/2.f;
//                 auto xpi = (x * M_PIf)/cutoff;
//                 b[k] = sin(xpi/factor)/xpi;
//             }
//             polyphase_[m] = FIRFilter(b);
//         }
//
//     }
//
//     int getFactor() const
//     {
//         if constexpr (FixedVar == PolyLength::Fixed) return M;
//         else
//             factor_;
//     }
//
//     template <class Ctxt, int N> void decimate(Ctxt c, delayline<N>
//     &delayline)
//     {
//         const auto id   = (getFactor() - c.getMultiRateId()) % getFactor();
//         auto &component = delayline.components_[id];
//         polyphase_[id].process(c, component[id]);
//
//         auto &x = c.getIn();
//         delayline.accumulator_ += x;
//
//         if (id == 0) {
//             x                = delayline.accumulator_;
//             delayline.accumulator_ = {0};
//         } else {
//             x = {0};
//         }
//     }
//
//     template <class Ctxt, int N> void interpolate(Ctxt c, delayline<N>
//     &delayline)
//     {
//         auto id = delayline.id_;
//         auto &x = c.getIn();
//
//         if (id == 0) {
//             delayline.accumulator_ = x;
//         } else {
//             x = delayline.accumulator_;
//         }
//         auto &component = delayline.components_[id];
//         polyphase_[id].process(c, component[id]);
//
//         ++delayline.id_;
//         if (delayline.id_ == getFactor()) {
//             delayline.id_ = 0;
//         }
//     }
//
//   private:
//     int id_{0};
//     std::enable_if_t<FixedVar == PolyLength::Var, int> factor_{M};
//     std::array<FIRFilter<Order>, M> polyphase_;
// };

} // namespace dsp
