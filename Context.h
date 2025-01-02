#pragma once

#include "Signal.h"

namespace dsp
{
template <int N, int BufferSize> class Context
{
  public:
    static constexpr auto Mask = BufferSize - 1;

    Context(Signal<N> *in, Signal<N> *buffer) : in_(in), buffer_(buffer) {}

    const Signal<N> &read(int i) const { return buffer_[(bufId_ + i) & Mask]; }
    void write(int i, const Signal<N> &x) { buffer_[(bufId_ + i) & Mask] = x; }

    Signal<N> &__restrict getIn() { return *in_; }

    void next(int incr = 1)
    {
        nextBufId(incr);
        nextIn(incr);
    }
    void nextBlock() { next(blockSize_); }

    void setIn(Signal<N> *in) { in_ = in; }
    void setBlockSize(int blockSize) { blockSize_ = blockSize; }
    int getBlockSize() const { return blockSize_; }

    constexpr int getIncr() const { return 1; }

  protected:
    void nextBufId(int incr) { bufId_ = (bufId_ - incr) & Mask; }
    void nextIn(int incr) { in_ += incr; }

  private:
    int bufId_{0};
    Signal<N> *buffer_;

    int blockSize_;
    Signal<N> *in_;
};

template <class Ctxt, typename P> void _processBlock(Ctxt c, P process)
{
    for (int n = 0; n < c.getBlockSize(); n += c.getIncr()) {
        process(c);
        c.next();
    }
}

// macro to easily redefine
#define processBlock(c, capt, func) _processBlock(c, capt(decltype(c) c) func)
} // namespace dsp
