#pragma once

#include <cassert>
#include <tuple>

namespace dsp
{
template <typename In, bool useVector = false> class Context
{
  public:
    Context(In *in, int blockSize) : blockSize_(blockSize), in_(in) {}
    Context(const Context &ctxt) = default;

    static constexpr auto VecSize       = useVector ? In::VectorSize : 1;
    static constexpr auto isUsingVector = useVector;
    using Type = typename std::conditional<useVector, typename In::Vector,
                                           typename In::Scalar>::type;

    auto vec() const { return Context<In, true>(in_, blockSize_); }
    int vecSize() const { return VecSize; }
    int getStep() const { return VecSize; }

    auto &__restrict getIn()
    {
        if constexpr (useVector) {
            return in_->toVector();
        } else {
            return in_->toScalar();
        }
    }
    void setIn(typename In::Scalar &x) { in_ = &x[0]; }
    void setIn(typename In::Vector &x) { in_ = &x[0]; }

    void next(int incr = VecSize) { nextIn(incr); }

    void nextBlock() { next(blockSize_); }
    int getBlockSize() const { return blockSize_; }

  protected:
    void nextIn(int incr) { in_ += incr; }

  private:
    const int blockSize_;
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

#define contextFor(ctxt)                                          \
    for (auto [c, n] = std::tuple{ctxt, 0}; n < c.getBlockSize(); \
         n += c.getStep(), c.next())

#define vecFor(ctxt) for (auto i = 0; i < ctxt.getStep(); ++i)

// macro to easily redefine
#define COMMA                 ,
#define processBlock(c, func) _processBlock([&](decltype(c) c) func, c)
#define processBlock2(c1, c2, func) \
    _processBlock([&](decltype(c1) c1, decltype(c2) c2) func, c1, c2)
} // namespace dsp
