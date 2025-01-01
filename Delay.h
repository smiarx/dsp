#pragma once

#include "Signal.h"
#include <cassert>

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

template <int N, int BuffSize> class Context
{
  public:
    static constexpr auto Mask = BuffSize - 1;

    Context(Signal<N> *buffer) : buffer_(buffer) {}

    const Signal<N> &read(int i) const { return buffer_[(id_ + i) & Mask]; }
    void write(int i, Signal<N> &x) { buffer_[(id_ + i) & Mask] = x; }
    void next() { id_ -= 1; }

  private:
    int id_{0};
    Signal<N> *buffer_;
};

template <int N, int Length> class DelayLine
{
  public:
    DelayLine() : offset_(bufferOffset.add(Length)) {}

    template <int BufferSize> void write(Context<N, BufferSize> c, Signal<N> x)
    {
        c.write(offset_, x);
    }

    template <int BufferSize> Signal<N> read(Context<N, BufferSize> c, int id)
    {
        return c.read(offset_ + id);
    }

    template <int BufferSize> Signal<N> tail(Context<N, BufferSize> c)
    {
        return read(c, Length);
    }

  private:
    const int offset_;
};

template <int N> class DelayLine<N, 1>
{
  public:
    template <int BufferSize> void write(Context<N, BufferSize> c, Signal<N> x)
    {
        mem_ = x;
    }

    template <int BufferSize> Signal<N> read(Context<N, BufferSize> c, int i)
    {
        return mem_;
    }

    template <int BufferSize> Signal<N> tail(Context<N, BufferSize> c)
    {
        return mem_;
    }

  private:
    Signal<4> mem_;
};

template <int N, int L = 1> struct TapFix {

    template <int BufferSize, int Ld>
    Signal<N> read(Context<N, BufferSize> c, DelayLine<N, Ld> &d)
    {
        static_assert(L <= Ld, "tap delay length is bigger than delay line");
        return d.tail(c);
    }
};
template <int N> class TapNoInterp
{
  public:
    void setDelay(int i, int _id) { id_[i] = _id; }

    template <int BufferSize, int L>
    Signal<N> read(Context<N, BufferSize> c, DelayLine<N, L> &d)
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
  public:
    void setDelay(int i, float d)
    {
        auto id  = static_cast<int>(d);
        float fd = d - static_cast<float>(id);
        fd_[i]   = fd;
        setDelay(id);
    };

    template <int BufferSize, int L>
    Signal<N> read(Context<N, BufferSize> c, DelayLine<N, L> &d)
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

template <int N> class TapCub : public TapLin<N>
{
};

/*
template <int N> class TapAllPass : public TapNoInterp<N>
{
  public:
    template <int _N, int L> friend class TapLine;

    static constexpr auto minfdelay = 0.1f;
    void set(int i, float _d)
    {
        auto _id  = static_cast<int>(_d - minfdelay);
        float _fd = _d - static_cast<float>(_id);
        a[i]      = (1 - _fd) / (1 + _fd);
        TapNoInterp<N>::set(i, _id);
    };

    void set(Signal<N> _d)
    {
        for (int i = 0; i < N; ++i) {
            set(i, _d[i]);
        }
    }

  private:
    Signal<N> a;
    Signal<N> s;
};

template <template <int N> class D, int L> struct DelayType {
    template <int N> using Type  = D<N>;
    static constexpr auto Length = L;
};
*/

} // namespace dsp
