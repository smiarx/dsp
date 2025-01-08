#pragma once

#include "Delay.h"
#include "Signal.h"
#include <cmath>

namespace dsp
{
template <int Order> class FIRFilter
{
  public:
    static constexpr auto NCoeff = Order + 1;
    template <int N> using Mem   = CopyDelayLine<N, Order>;

    FIRFilter() = default;
    FIRFilter(std::array<float, NCoeff> &&b) : b_(b) {}

    template <class Ctxt, class Mem> void process(Ctxt c, Mem &mem) const
    {
        auto &x = c.getIn();

        auto out = x * b_[0];
        for (int n = 1; n < NCoeff; ++n) {
            auto xn = mem.read(c, n);
            out += xn * b_[n];
        }

        mem.write(c, x);
        x = out;
    }

  private:
    std::array<float, NCoeff> b_ = {0};
};

enum class PolyLength {
    Fixed,
    Var,
};

template <int Order, int M, PolyLength FixedVar = PolyLength::Fixed>
class FIRPolyphase
{
  public:
    template <class Ctxt> class MultiRateContext : public Ctxt
    {
        friend FIRPolyphase;

        MultiRateContext(Ctxt c, FIRPolyphase &ply) :
            Ctxt(c), ply_(ply), multiRateId_(ply_.id)
        {
        }

        void next(int incr = 1)
        {
            if (multiRateId_ == ply_.getFactor()) {
                auto incrNewRate = multiRateId_ / getFactor();
                Ctxt::next(incrNewRate);
                Ctxt::nextIn(incr - incrNewRate);

                multiRateId_ %= ply_.getFactor();
            } else {
                Ctxt::nextIn(incr);
            }
        }

        int getMultiRateId() const { return multiRateId_; }

        void save() { ply_.id_ = multiRateId_; }

      private:
        const FIRPolyphase &ply_;
        int multiRateId_;
    };

    template <class Ctxt>
    MultiRateContext<Ctxt> getMultiRateContext(Ctxt c) const
    {
        return MultiRateContext(c);
    }

    template <class Ctxt>
    Ctxt getContext(Ctxt c) { return Ctxt(c, getFactor()); }

    template <int N> struct Mem {
        Signal<N> accumulator_;
        std::array<typename FIRFilter<Order>::template Mem<N>, M> components_;
    };

    /* construct polyphase filter */
    FIRPolyphase(int factor, float cutoff = 1.f)
    {
        int N = factor*Order;
        for(int m = 0; m < M; ++m)
        {
            Signal<Order> b;
            for(int k = 0; k < Order; ++k)
            {
                int n = m+k*M;
                auto x = n - N/2.f;
                auto xpi = (x * M_PIf)/cutoff;
                b[k] = sin(xpi/factor)/xpi;
            }
            polyphase_[m] = FIRFilter(b);
        }
        
    }

    int getFactor() const
    {
        if constexpr (FixedVar == PolyLength::Fixed) return M;
        else
            factor_;
    }

    template <class Ctxt, int N> void decimate(Ctxt c, Mem<N> &mem)
    {
        const auto id   = (getFactor() - c.getMultiRateId()) % getFactor();
        auto &component = mem.components_[id];
        polyphase_[id].process(c, component[id]);

        auto &x = c.getIn();
        mem.accumulator_ += x;

        if (id == 0) {
            x                = mem.accumulator_;
            mem.accumulator_ = {0};
        } else {
            x = {0};
        }
    }

    template <class Ctxt, int N> void interpolate(Ctxt c, Mem<N> &mem)
    {
        auto id = mem.id_;
        auto &x = c.getIn();

        if (id == 0) {
            mem.accumulator_ = x;
        } else {
            x = mem.accumulator_;
        }
        auto &component = mem.components_[id];
        polyphase_[id].process(c, component[id]);

        ++mem.id_;
        if (mem.id_ == getFactor()) {
            mem.id_ = 0;
        }
    }

  private:
    int id_{0};
    std::enable_if_t<FixedVar == PolyLength::Var, int> factor_{M};
    std::array<FIRFilter<Order>, M> polyphase_;
};

} // namespace dsp
