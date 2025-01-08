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

template <int L = 0, int Off = 0> class DelayLine;

template <> class DelayLine<>
{
  public:
    /*
     * Implement a delay line
     * templates:
     *      N: vector size
     *      Length: maximum length of the delay line
     */

    static constexpr int Length = 10000000;
    static constexpr int NextOffset = 0;

    template<int>
    using WithOffset = DelayLine;

    DelayLine(int length) : offset_(bufferOffset.add(length)) {}

    /* write value inside delay line */
    template <class Ctxt, class In> void write(Ctxt c, const In &x) const
    {
        c.getBuffer().write(offset_, x);
    }

    /* read value at delay id */
    template <class Ctxt> const auto &read(Ctxt c, int id) const
    {
        return c.getBuffer().read(offset_ + id);
    }

    template <class Ctxt, int Size> auto readContiguous(Ctxt c, int id) const
    {
        return c.getBuffer().template readContiguous<Size>(offset_ + id);
    }

  private:
    const int offset_;
};

/* delayline with compile time length */
template <int L, int Off> class DelayLine
{
  public:
    static constexpr auto Length = L;
    static constexpr auto Offset = Off;
    static constexpr auto NextOffset = Offset+Length;

    template<int O>
    using WithOffset = DelayLine<L,O>;

    constexpr DelayLine() {}

    /* write value inside delay line */
    template <class Ctxt, class In> void write(Ctxt c, const In &x) const
    {
        c.getBuffer().write(Offset, x);
    }

    /* read value at delay id */
    template <class Ctxt> const auto &read(Ctxt c, int id) const
    {
        return c.getBuffer().read(Offset + id);
    }

    /* read value at tail */
    template <class Ctxt> const auto &tail(Ctxt c) const
    {
        return read(c, Length);
    }
};

template <int N, int L = 1, int Off=0> class CopyDelayLine
{
    /* */
  public:
    static constexpr auto Length = L;
    static constexpr auto Offset = Off;
    static constexpr auto NextOffset = Offset;

    template<int O>
    using WithOffset = CopyDelayLine<N,L,O>;

    template <class Ctxt> void write(Ctxt c, const Signal<N> &x)
    {
        for (int j = Length-1; j > 0; --j) mem_[j] = mem_[j-1];
        mem_[0] = x;
    }

    template <class Ctxt> const Signal<N> &read(Ctxt c, int i) const
    {
        return mem_[i - 1];
    }

    template <class Ctxt, int Size> auto readContiguous(Ctxt c, int i) const
    {
        std::array<Signal<N>, Size> val = mem_+i;
        return val;
    }

    template <class Ctxt> const Signal<N> &tail(Ctxt c) const
    {
        return mem_[Length - 1];
    }

  private:
    Signal<N> mem_[Length] = {0};
};

template <class DL, class DLi, int Off = 0> class NestedDelayLine : public DL::template WithOffset<Off>
{
    using Outer = typename DL::template WithOffset<Off>;
    using Inner = typename DLi::template WithOffset<Outer::NextOffset>;
  public:

    static constexpr auto Length = Outer::Length+Inner::Length;
    static constexpr auto Offset = Off;
    static constexpr auto NextOffset = Inner::NextOffset;

    template<int O>
    using WithOffset = NestedDelayLine<DL,DLi,O>;

    NestedDelayLine() = default;
    NestedDelayLine(Outer &&dl, Inner &&dli) : Outer(dl), inner_(dli) {}

    Inner& getInner() { return inner_;}

  private:
    Inner inner_;
};

template <class DL, int Nm, int Off = 0> class ArrayDelayLine :
    protected std::array<typename DL::template WithOffset<Off>, Nm>
{
    using Base = typename DL::template WithOffset<Off>;
public:
    static constexpr auto Length = DL::Length*Nm;
    static constexpr auto Offset = Off;
    static constexpr auto OneOffset = (Base::NextOffset-Offset);
    static constexpr auto NextOffset = OneOffset*Nm;

    template<int O>
    using WithOffset = ArrayDelayLine<DL,Nm,O>;

    template <class Ctxt>
    auto& get(Ctxt& c, int i)
    {
        c.nextBufId(OneOffset*i);
        return (*this)[i];
    }
};

#define nextTo(delayline) decltype(delayline)::NextOffset

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

/* help function */
template <class Ctxt, class DL, class T>
void _fixread(Ctxt c, const DL &delayline, T& x, int i) {}
template <class Ctxt, class DL, class T, int D, int...Ds>
void _fixread(Ctxt c, const DL &delayline, T& x, int i)
{
    static_assert(D <= DL::Length,
                  "tap delay length is bigger than delay line");
    x[i] = delayline.read(c, D)[i];
    _fixread<Ctxt, DL, T, Ds...>(c, delayline, x, ++i);
}
template <int D=1, int...Ds> struct TapFix : public TapFix<Ds...>{
    template <class Ctxt, class DL>
    auto read(Ctxt c, const DL &delayline) const
    {
        std::remove_const_t<std::remove_reference_t<decltype(delayline.read(c,0))>> x = {0};
        _fixread<Ctxt, DL, decltype(x), D, Ds...>(c, delayline, x, 0);
        return x;
    }
};
template <int D> struct TapFix<D> {
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
