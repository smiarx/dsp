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
    Context(const Context &ctxt)            = default;
    Context &operator=(const Context &ctxt) = default;

    static constexpr auto VecSize       = Vectorize ? In::VectorSize : 1;
    static constexpr auto isUsingVector = Vectorize;
    using BaseType                      = In;
    using Type =
        std::conditional_t<Vectorize, typename In::Vector, typename In::Scalar>;

    [[nodiscard]] auto vec() const
    {
        return Context<In, true>(in_, blockSize_);
    }
    [[nodiscard]] auto scalar() const
    {
        return Context<In, false>(in_, blockSize_);
    }

    [[nodiscard]] int vecSize() const { return VecSize; }
    [[nodiscard]] int getStep() const { return VecSize; }

    auto &getSignal() { return in_->template toSignal<Vectorize>(); }

    template <typename T> void setSamples(T &x) { in_ = &x[0]; }

    void next(int incr = VecSize) { nextIn(incr); }

    [[nodiscard]] int getBlockSize() const { return blockSize_; }
    void setBlockSize(int blockSize) { blockSize_ = blockSize; }

  protected:
    void nextIn(int incr) { in_ += incr; }

  private:
    int blockSize_;
    In *__restrict in_;
};

template <class Ctxt, class... Ctxts>
void _ctxtInfos(int &blockSize, int &step, Ctxt c, Ctxts... cs)
{
    blockSize = c.getBlockSize();
    step      = c.getStep();
    /* check if blockSize and incr are the same */
    (
        [&] {
            auto _blockSize = cs.getBlockSize();
            assert(blockSize == _blockSize);
        }(),
        ...);
    (
        [&] {
            auto _step = cs.getStep();
            assert(step == _step);
        }(),
        ...);
}

#define contextFor(ctxt)                                          \
    for (auto [c, n] = std::tuple{ctxt, 0}; n < c.getBlockSize(); \
         n += c.getStep(), c.next())

#define vecFor(ctxt) for (auto i = 0; i < ctxt.getStep(); ++i)

#define arrayFor(x, k) for (size_t k = 0; k < x.size(); ++k)
#define inFor(x, k, i)                    \
    for (size_t k = 0; k < x.size(); ++k) \
        for (size_t i = 0; i < x[0].size(); ++i)

#define PROCESSBLOCK_                                \
    template <class Ctxt, class State>               \
    void processBlock(Ctxt ctxt, State &state) const \
    {                                                \
        contextFor(ctxt) { process(c, state); }      \
    }

} // namespace dsp
