#pragma once

#include <array>
#include <cassert>

namespace dsp
{
template <typename In> class Context
{
  public:
    Context(In *in) : in_(in) {}

    In &__restrict getIn() { return *in_; }

    void next(int incr = 1) { nextIn(incr); }

    void nextBlock() { next(blockSize_); }

    void setIn(In *in) { in_ = in; }
    void setBlockSize(int blockSize) { blockSize_ = blockSize; }
    int getBlockSize() const { return blockSize_; }
    constexpr int getIncr() const { return 1; }

  protected:
    void nextIn(int incr) { in_ += incr; }
    void nextBufId(int incr) {}

  private:
    int blockSize_;
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
#define processBlock(c, capt, func) _processBlock(capt(decltype(c) c) func, c)
#define processBlock2(c1, c2, capt, func) \
    _processBlock(capt(decltype(c1) c1, decltype(c2) c2) func, c1, c2)
} // namespace dsp
