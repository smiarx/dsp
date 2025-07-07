#pragma once

#include "FIRFilter.h"

namespace dsp
{

template <typename T, size_t Order, size_t M>
class MultiRate : public MultiRate<T, Order, M - 1>
{
  public:
    using Base = MultiRate<T, Order, 1>;
    using Next = MultiRate<T, Order, M - 1>;

    template <size_t Offset = 0>
    using DLDecimate = typename FIRDecimate<T, Order, M>::template DL<Offset>;
    template <size_t Offset = 0>
    using DLInterpolate =
        typename FIRInterpolate<T, Order, M>::template DL<Offset>;

    MultiRate(baseType<T> cutoff = 1) :
        firdecimate_(cutoff), firinterpolate_(cutoff)
    {
    }

    template <class CtxtIn, class CtxtOut, class DL>
    int decimate(unsigned int rate, CtxtIn cin, CtxtOut &cout, DL &dl,
                 int decimateId) const
    {
        if (rate == M)
            return firdecimate_.decimate(cin.vec(), cout, dl, decimateId);
        else
            return Next::decimate(rate, cin, cout, dl, decimateId);
    }

    template <class CtxtIn, class CtxtOut, class DL>
    int interpolate(unsigned int rate, CtxtIn cin, CtxtOut &cout, DL &dl,
                    int interpolateId) const
    {
        if (rate == M)
            return firinterpolate_.interpolate(cin, cout, dl, interpolateId);
        else
            return Next::interpolate(rate, cin, cout, dl, interpolateId);
    }

  private:
    const FIRDecimate<T, Order, M> firdecimate_;
    const FIRInterpolate<T, Order, M> firinterpolate_;
};

template <typename T, size_t Order> class MultiRate<T, Order, 1>
{
  public:
// function seems more efficient with unknown behaviour
// when rate > M. be careful
// ignore warning "return-type'
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#elif defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4715)
#endif

    template <typename F, bool Vec, class DL>
    int decimate(unsigned int rate, Context<F, Vec> cin, Context<F, Vec> &cout,
                 DL &, int) const
    {
        if (rate == 1) {
            cout = cin;
            return 0;
        } else
            assert(false);
    }
    template <typename F, bool Vec, class DL>
    int interpolate(unsigned int rate, Context<F, Vec> cin,
                    Context<F, Vec> &cout, DL &, int) const
    {
        if (rate == 1) {
            cout = cin;
            return 0;
        } else
            assert(false);
    }

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#pragma warning(pop)
#endif

    constexpr int getDelay(unsigned int rate) { return (Order + 1) * rate - 1; }
};

} // namespace dsp
