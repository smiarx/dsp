#pragma once

#include "Filter.h"
#include "Signal.h"

namespace dsp
{
/* scipy.signal.cheby1(10,2,1,analog=True,output='sos') */
constexpr float decimateFilterBA[][2][3] = {
    {{0., 0., 0.00255383}, {1., 0.21436212, 0.0362477}},
    {{0., 0., 1.}, {1., 0.19337886, 0.21788333}},
    {{0., 0., 1.}, {1., 0.15346633, 0.51177596}},
    {{0., 0., 1.}, {1., 0.09853145, 0.80566858}},
    {{0., 0., 1.}, {1., 0.03395162, 0.98730422}}};

template <int N, int BufferSize>
class DownSampleContext : public Context<N, BufferSize>
{
    using Ctxt = Context<N, BufferSize>;

  public:
    DownSampleContext(Ctxt c, int decimate, int decimateId) :
        Ctxt(c), decimate_(decimate)
    {
        Ctxt::nextIn(decimateId);
    }

    int getIncr() const { return decimate_; }

    void next(int incr = 1)
    {
        Ctxt::nextBufId(incr);
        Ctxt::nextIn(incr * decimate_);
    }

  private:
    int decimate_;
    int decimateId_;
};

template <int N> class DownSampler
{
  public:
    template <class Ctxt> void process(Ctxt c)
    {
        auto &x = c.getIn();

        if (decimate_ == 1) return;

        lpFilter.process(c, decimateFilterDL);

        if (decimateId_ == 0) {
            for (int i = 0; i < N; ++i) x[i] *= decimate_;
        } else {
            x = 0;
        }
        decimateId_ = (decimateId_ + 1) % decimate_;
    }

    template <class Ctxt> auto getContext(Ctxt c)
    {
        return DownSampleContext(c, decimate_, decimateId_);
    }

    void set(int decimate)
    {
        decimate_ = decimate;

        if (decimate_ == 1) {
            decimateId_ = 0;
        } else {
            Signal<1> f{1.f / decimate_};
            lpFilter.sosanalog(decimateFilterBA, f);
        }
    }

  private:
    int decimate_{1};
    int decimateId_{0};
    Filter<N, sizeof(decimateFilterBA) / sizeof(decimateFilterBA[0])> lpFilter;
    typename decltype(lpFilter)::DelayLineType decimateFilterDL;
};

} // namespace dsp
