#pragma once

#include <cassert>

namespace dsp
{
template <typename In> class Context
{
  public:
    Context(Context &) = default;
    Context(In *in, int blockSize, int step = 1) :
        blockSize_(blockSize), step_(step), in_(in)
    {
    }
    Context(Context &ctxt, int step = 1) :
        Context(ctxt), blockSize_(blockSize_ / step), step_(step)
    {
    }

    In &__restrict getIn() { return *in_; }

    void next(int incr = 1) { nextIn(incr); }

    void nextBlock() { next(blockSize_); }

    void setIn(In *in) { in_ = in; }
    int getBlockSize() const { return blockSize_; }
    int getStep() const { return step_; }

  protected:
    void nextIn(int incr) { in_ += incr; }

  private:
    const int blockSize_;
    const int step_{1};
    In *in_;
};

template <class Ctxt, class... Ctxts>
void _ctxtInfos(int &blockSize, int &incr, Ctxt c, Ctxts... cs)
{
    blockSize = c.getBlockSize();
    incr      = c.getIncr();
    /* check if blockSize and incr are the same */
    (
        [&] {
            auto _blockSize = cs.getBlockSize();
            assert(blockSize == _blockSize);
        }(),
        ...);
    (
        [&] {
            auto _incr = cs.getIncr();
            assert(incr == _incr);
        }(),
        ...);
}
template <typename P, class... Ctxt> void _processBlock(P process, Ctxt... cs)
{
    int blockSize, incr;
    _ctxtInfos(blockSize, incr, cs...);
    for (int n = 0; n < blockSize; n += incr) {
        process(cs...);
        (
            [&] {
                cs.next();
            }(),
            ...);
    }
}

// macro to easily redefine
#define COMMA                 ,
#define processBlock(c, func) _processBlock([&](decltype(c) c) func, c)
#define processBlock2(c1, c2, func) \
    _processBlock([&](decltype(c1) c1, decltype(c2) c2) func, c1, c2)
} // namespace dsp
