#pragma once

#include <cassert>

#include "MultiVal.h"

namespace dsp
{
template <typename T, bool Vec = false> class Context
{
  public:
    Context(T *data, int blockSize = 1) : blockSize_(blockSize), data_(data) {}
    Context(const Context &ctxt)            = default;
    Context &operator=(const Context &ctxt) = default;

    static constexpr auto kIncrSize = Vec ? sizeof(batch<T>) / sizeof(T) : 1;
    static constexpr bool kUseVec   = Vec;

    using BaseType = T;
    using Type     = std::conditional_t<Vec, batch<T>, T>;
    // input and output type
    using SigType = std::conditional_t<Vec, decltype(loadfuncs::loadBatch(T{})),
                                       decltype(loadfuncs::load(T{}))>;

    [[nodiscard]] auto load(const T &x) const
    {
        if constexpr (Vec) {
            return loadfuncs::loadBatch(x);
        } else {
            return loadfuncs::load(x);
        }
    }

    [[nodiscard]] auto store(T &dest, SigType val) const
    {
        if constexpr (Vec) {
            return loadfuncs::storeBatch(dest, val);
        } else {
            return loadfuncs::store(dest, val);
        }
    }

    [[nodiscard]] auto getInput() const { return load(*data_); }
    template <typename S> void setOutput(S value) { store(*data_, value); }

    void next(int incr = kIncrSize) { nextData(incr); }

    [[nodiscard]] int getBlockSize() const { return blockSize_; }
    void setBlockSize(int blockSize) { blockSize_ = blockSize; }

    void setData(T *data) { data_ = data; }
    [[nodiscard]] const T *getData() const { return data_; }

    [[nodiscard]] auto vec() const
    {
        return Context<T, true>(data_, blockSize_);
    }
    [[nodiscard]] auto scalar() const
    {
        return Context<T, false>(data_, blockSize_);
    }

  protected:
    void nextData(int incr) { data_ += incr; }

  private:
    int blockSize_;
    T *__restrict data_;
};

template <class Ctxt> class ContextRun
{
  public:
    ContextRun(Ctxt ctxt) : ctxt_(std::move(ctxt)) {}

    template <typename Func>
    ContextRun(Ctxt ctxt, Func func) : ContextRun(std::move(ctxt))
    {
        run(func);
    }

    template <typename Func> void run(Func func)
    {
        auto n                  = 0;
        constexpr int kIncrSize = Ctxt::kIncrSize;
        for (; n < ctxt_.getBlockSize() - kIncrSize + 1;
             n += kIncrSize, ctxt_.next()) {
            func(ctxt_);
        }

        if constexpr (Ctxt::kUseVec) {
            auto ctxtScal = ctxt_.scalar();

            for (; n < ctxtScal.getBlockSize(); n += 1, ctxtScal.next()) {
                func(ctxtScal);
            }
        }
    }

    template <typename Func> ContextRun &operator=(Func func)
    {
        run(func);
        return *this;
    }

    explicit operator bool() { return true; }

  private:
    Ctxt ctxt_;
};

#define CTXTRUN(ctxt) \
    if (dsp::ContextRun contextRun{ctxt}) contextRun = [&](auto ctxt)
#define CTXTRUNVEC(ctxt)                           \
    if (dsp::ContextRun contextRunVec{ctxt.vec()}) \
    contextRunVec = [&](auto ctxt)

#define PROCESSBLOCK_                                \
    template <class Ctxt, class State>               \
    void processBlock(Ctxt ctxt, State &state) const \
    {                                                \
        CTXTRUN(ctxt) { process(ctxt, state); };     \
    }

using namespace loadfuncs;

} // namespace dsp
