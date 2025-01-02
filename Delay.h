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

template <int N, int Length, class V = void> class DelayLine
{
    /*
     * Implement a delay line
     * templates:
     *      N: vector size
     *      Length: maximum length of the delay line
     */
  public:
    DelayLine() : offset_(bufferOffset.add(Length)) {}

    /* write value inside delay line */
    template <class Ctxt> void write(Ctxt c, const Signal<N> &x)
    {
        c.getBuffer().write(offset_, x);
    }

    /* read value at delay id */
    template <class Ctxt> const Signal<N> &read(Ctxt c, int id) const
    {
        return c.getBuffer().read(offset_ + id);
    }

    /* read value at tail */
    template <class Ctxt> const Signal<N> &tail(Ctxt c) const
    {
        return read(c, Length);
    }

  private:
    const int offset_;
};

constexpr auto maxCopyDelay = 4;
template <int N, int Length>
class DelayLine<N, Length, std::enable_if_t<Length <= maxCopyDelay>>
{
    /* if Length is <= 4 we implement the delayline with a copy system instead
     * of a ringbuffer */
  public:
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

/*
 * Taps: different classes that can read into a delay line with different
 * position and interpolation In each class DelayLineType<L> define the type of
 * the base DelayLine that can be read by the tap
 */

template <int N> struct TapTail {
    /* Tap that reads a the end of a delayline */
    template <int L> using DelayLineType = DelayLine<N, L>;

    template <class Ctxt, int L>
    const Signal<N> &read(Ctxt c, DelayLine<N, L> &d) const
    {
        return d.tail(c);
    }
};

template <int N, int D = 1> struct TapFix {
    /* Tap that reads at a fix point (delay D) in the delay line */
    template <int...> using DelayLineType = DelayLine<N, D>;

    template <class Ctxt, int L>
    const Signal<N> &read(Ctxt c, const DelayLine<N, L> &d) const
    {
        static_assert(D <= L, "tap delay length is bigger than delay line");
        return d.read(c, D);
    }
};
template <int N> class TapNoInterp
{
    /* Tap that can read at any point in de delay line, without interpolation */
  public:
    template <int... D> using DelayLineType = DelayLine<N, D...>;

    void setDelay(int i, int id) { id_[i] = id; }
    void setDelay(iSignal<N> id) { id_ = id; }

    template <class Ctxt, int L>
    Signal<N> read(Ctxt c, DelayLine<N, L> &d) const
    {
        Signal<N> x;
        for (int i = 0; i < N; i++) {
            assert(id_[i] <= L);
            x[i] = d.read(c, id_[i])[i];
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
    template <int D> using DelayLineType = DelayLine<N, D>;

    void setDelay(int i, float d)
    {
        auto id  = static_cast<int>(d);
        float fd = d - static_cast<float>(id);
        fd_[i]   = fd;
        setDelay(id);
    };

    template <class Ctxt, int L>
    Signal<N> read(Ctxt c, DelayLine<N, L> &d) const
    {
        Signal<N> x0;
        Signal<N> x1;

        auto id = TapNoInterp<N>::id_;

        for (int i = 0; i < N; ++i) {
            assert(id[i] < L);
            x0[i] = d.read(c, id[i])[i];
            x1[i] = d.read(c, id[i] + 1)[i];
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
