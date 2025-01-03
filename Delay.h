#pragma once

#include "Context.h"
#include "Signal.h"
#include <cassert>
#include <type_traits>

namespace dsp
{

class BufferOffset
{
  public:
    int add(int length)
    {
        const auto offset = value_;
        value_ += length + 1;
        return offset;
    }

  private:
    int value_{0};
};

static BufferOffset bufferOffset;

template <int L = 0> class DelayLine;

template <> class DelayLine<0>
{
  public:
    static constexpr unsigned int Length = -1;
    /*
     * Implement a delay line
     * templates:
     *      N: vector size
     *      Length: maximum length of the delay line
     */
    DelayLine(int length) : offset_(bufferOffset.add(length)) {}

    /* write value inside delay line */
    template <class Ctxt, class In> void write(Ctxt c, const In &x)
    {
        c.getBuffer().write(offset_, x);
    }

    /* read value at delay id */
    template <class Ctxt> const auto &read(Ctxt c, int id) const
    {
        return c.getBuffer().read(offset_ + id);
    }

  private:
    const int offset_;
};

/* delayline with compile time length */
template <int L> class DelayLine : public DelayLine<0>
{
  public:
    static constexpr auto Length = L;

    DelayLine() : DelayLine<0>(Length) {}

    /* read value at tail */
    template <class Ctxt> const auto &tail(Ctxt c) const
    {
        return read(c, Length);
    }
};

template <int N, int L = 1> class CopyDelayLine
{
    /* */
  public:
    static constexpr auto Length = L;

    template <class Ctxt> void write(Ctxt c, const Signal<N> &x)
    {
        for (int j = 0; j < Length - 1; ++j) mem_[j + 1] = mem_[j];
        mem_[0] = x;
    }

    template <class Ctxt> const Signal<N> &read(Ctxt c, int i) const
    {
        return mem_[i - 1];
    }

    template <class Ctxt> const Signal<N> &tail(Ctxt c) const
    {
        return mem_[Length - 1];
    }

  private:
    Signal<N> mem_[Length] = {0};
};

template <class DL, class DLi> class NestedDelayLine : public DL
{
  public:
    NestedDelayLine() = default;
    NestedDelayLine(DL &&dl, DLi &&dli) : DL(dl), inner_(dli) {}
    DLi inner_;
};

/*
 * Taps: different classes that can read into a delay line with different
 * position and interpolation In each class DelayLineType<L> define the type of
 * the base DelayLine that can be read by the tap
 */

struct TapTail {
    /* Tap that reads at the end of a delayline */

    template <class Ctxt, class DL>
    const auto &read(Ctxt c, const DL &delayline) const
    {
        return delayline.tail(c);
    }
};

template <int D = 1> struct TapFix {
    /* Tap that reads at a fix point (delay D) in the delay line */

    template <class Ctxt, class DL>
    const auto &read(Ctxt c, const DL &delayline) const
    {
        static_assert(D <= DL::Length,
                      "tap delay length is bigger than delay line");
        return delayline.read(c, D);
    }
};
template <int N> class TapNoInterp
{
    /* Tap that can read at any point in de delay line, without interpolation */
  public:
    void setDelay(int i, int id) { id_[i] = id; }
    void setDelay(iSignal<N> id) { id_ = id; }

    template <class Ctxt, class DL>
    Signal<N> read(Ctxt c, const DL &delayline) const
    {
        Signal<N> x;
        for (int i = 0; i < N; i++) {
            assert(id_[i] <= DL::Length);
            x[i] = delayline.read(c, id_[i])[i];
        }
        return x;
    }

  protected:
    iSignal<N> id_;
};

template <int N> class TapLin : public TapNoInterp<N>
{
    /* Tap that can read at any point in de delay line, with linear
     * interpolation */
  public:
    void setDelay(int i, float d)
    {
        auto id  = static_cast<int>(d);
        float fd = d - static_cast<float>(id);
        fd_[i]   = fd;
        setDelay(id);
    };

    template <class Ctxt, class DL> Signal<N> read(Ctxt c, DL &delayline) const
    {
        Signal<N> x0;
        Signal<N> x1;

        auto id = TapNoInterp<N>::id_;

        for (int i = 0; i < N; ++i) {
            assert(id[i] < DL::Length);
            x0[i] = delayline.read(c, id[i])[i];
            x1[i] = delayline.read(c, id[i] + 1)[i];
        }
        for (int i = 0; i < N; ++i) {
            x0[i] += fd_[i] * (x1[i] - x0[i]);
        }
        return x0;
    }

  protected:
    Signal<N> fd_;
};
} // namespace dsp
