#pragma once

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

template <int N, int BuffSize> class Context
{
  public:
    static constexpr auto Mask = BuffSize - 1;

    Context(Signal<N> *buffer) : buffer_(buffer) {}

    const Signal<N> &read(int i) const { return buffer_[(id_ + i) & Mask]; }
    void write(int i, const Signal<N> &x) { buffer_[(id_ + i) & Mask] = x; }
    void next() { id_ -= 1; }

  private:
    int id_{0};
    Signal<N> *buffer_;
};

template <int N, int Length, class V = void> class DelayLine
{
  public:
    DelayLine() : offset_(bufferOffset.add(Length)) {}

    template <int BufferSize>
    void write(Context<N, BufferSize> c, const Signal<N> &x)
    {
        c.write(offset_, x);
    }

    template <int BufferSize>
    const Signal<N> &read(Context<N, BufferSize> c, int id) const
    {
        return c.read(offset_ + id);
    }

    template <int BufferSize>
    const Signal<N> &tail(Context<N, BufferSize> c) const
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
  public:
    template <int BufferSize>
    void write(Context<N, BufferSize> c, const Signal<N> &x)
    {
        for (int j = 0; j < Length; ++j) mem_[j + 1] = mem_[j];
        mem_[0] = x;
    }

    template <int BufferSize>
    const Signal<N> &read(Context<N, BufferSize> c, int i) const
    {
        return mem_[i - 1];
    }

    template <int BufferSize>
    const Signal<N> &tail(Context<N, BufferSize> c) const
    {
        return mem_[Length - 1];
    }

  private:
    Signal<N> mem_[Length];
};

template <int N> struct TapTail {
    template <int BufferSize, int L>
    const Signal<N> &read(Context<N, BufferSize> c, DelayLine<N, L> &d) const
    {
        return d.tail(c);
    }
};

template <int N, int D = 1> struct TapFix {

    template <int BufferSize, int L>
    const Signal<N> &read(Context<N, BufferSize> c,
                          const DelayLine<N, L> &d) const
    {
        static_assert(D <= L, "tap delay length is bigger than delay line");
        return d.read(c, D);
    }
};
template <int N> class TapNoInterp
{
  public:
    void setDelay(int i, int _id) { id_[i] = _id; }

    template <int BufferSize, int L>
    Signal<N> read(Context<N, BufferSize> c, DelayLine<N, L> &d) const
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
    Signal<N> read(Context<N, BufferSize> c, DelayLine<N, L> &d) const
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
} // namespace dsp
