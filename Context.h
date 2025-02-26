#pragma once

#include <cassert>
#include <cstdlib>
#include <tuple>

namespace dsp
{
template <typename In, bool Vectorize = false> class Context
{
  public:
    Context(In *in, int blockSize = 1) : blockSize_(blockSize), in_(in) {}
    Context(const Context &ctxt) = default;

    static constexpr auto VecSize       = Vectorize ? In::VectorSize : 1;
    static constexpr auto isUsingVector = Vectorize;
    using BaseType                      = In;
    using Type = typename std::conditional<Vectorize, typename In::Vector,
                                           typename In::Scalar>::type;

    auto vec() const { return Context<In, true>(in_, blockSize_); }
    auto scalar() const { return Context<In, false>(in_, blockSize_); }

    int vecSize() const { return VecSize; }
    int getStep() const { return VecSize; }

    auto &__restrict getSignal() { return in_->template toSignal<Vectorize>(); }

    template <typename T> void setSamples(T &x) { in_ = &x[0]; }

    void next(int incr = VecSize) { nextIn(incr); }

    void nextBlock() { next(blockSize_); }
    int getBlockSize() const { return blockSize_; }
    void setBlockSize(int blockSize) { blockSize_ = blockSize; }

  protected:
    void nextIn(int incr) { in_ += incr; }

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

#define contextFor(ctxt)                                          \
    for (auto [c, n] = std::tuple{ctxt, 0}; n < c.getBlockSize(); \
         n += c.getStep(), c.next())

#define vecFor(ctxt) for (auto i = 0; i < ctxt.getStep(); ++i)

#define arrayFor(x, k) for (size_t k = 0; k < x.size(); ++k)
#define inFor(x, k, i)                    \
    for (size_t k = 0; k < x.size(); ++k) \
        for (size_t i = 0; i < x[0].size(); ++i)

// macro to easily redefine
#define COMMA                 ,
#define processBlock(c, func) _processBlock([&](decltype(c) c) func, c)
#define processBlock2(c1, c2, func) \
    _processBlock([&](decltype(c1) c1, decltype(c2) c2) func, c1, c2)
} // namespace dsp
