#pragma once

#include <cassert>

#include "simd/multi.h"

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
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
    using SigType = std::conditional_t<Vec, decltype(dsp::loadBatch(T{})),
                                       decltype(dsp::load(T{}))>;

    void load(T &&) const { static_assert(false, "do not load rvalue"); }
    [[nodiscard]] auto load(const T &x) const
    {
        if constexpr (Vec) {
            return dsp::DSP_ARCH_NAMESPACE::loadBatch(x);
        } else {
            return dsp::DSP_ARCH_NAMESPACE::load(x);
        }
    }

    [[nodiscard]] auto store(T &dest, SigType val) const
    {
        if constexpr (Vec) {
            return dsp::DSP_ARCH_NAMESPACE::storeBatch(dest, val);
        } else {
            return dsp::DSP_ARCH_NAMESPACE::store(dest, val);
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

template <bool Rev, class Ctxt1, class Ctxt2> class ContextRun2
{
  public:
    ContextRun2(Ctxt1 ctxt1, Ctxt2 ctxt2) :
        ctxt1_(std::move(ctxt1)), ctxt2_(std::move(ctxt2))
    {
    }

    template <typename Func>
    ContextRun2(Ctxt1 ctxt1, Ctxt2 ctxt2, Func func) :
        ContextRun2(std::move(ctxt1), std::move(ctxt2))
    {
        run(func);
    }

    template <typename Func> void run(Func func)
    {
        static_assert(Ctxt1::kUseVec == Ctxt2::kUseVec);
        assert(ctxt1_.getBlockSize() == ctxt2_.getBlockSize());

        auto n                   = 0;
        constexpr int kIncrSize  = Ctxt1::kIncrSize;
        constexpr auto kNextStep = Rev ? -kIncrSize : kIncrSize;

        if constexpr (Rev) {
            ctxt1_.next(ctxt1_.getBlockSize() - kIncrSize);
            ctxt2_.next(ctxt2_.getBlockSize() - kIncrSize);
        }

        for (; n < ctxt1_.getBlockSize() - kIncrSize + 1;
             n += kIncrSize, ctxt1_.next(kNextStep), ctxt2_.next(kNextStep)) {
            func(ctxt1_, ctxt2_);
        }

        if constexpr (Ctxt1::kUseVec) {
            auto ctxt1Scal = ctxt1_.scalar();
            auto ctxt2Scal = ctxt2_.scalar();

            constexpr auto kNextStep = Rev ? -1 : 1;

            for (; n < ctxt1Scal.getBlockSize();
                 n += 1, ctxt1Scal.next(kNextStep), ctxt2Scal.next(kNextStep)) {
                func(ctxt1Scal, ctxt2Scal);
            }
        }
    }

    template <typename Func> ContextRun2 &operator=(Func func)
    {
        run(func);
        return *this;
    }

    explicit operator bool() { return true; }

  private:
    Ctxt1 ctxt1_;
    Ctxt2 ctxt2_;
};

#define CTXTRUN(ctxt) \
    if (dsp::ContextRun contextRun{ctxt}) contextRun = [&](auto ctxt)
#define CTXTRUNVEC(ctxt)                           \
    if (dsp::ContextRun contextRunVec{ctxt.vec()}) \
    contextRunVec = [&](auto ctxt)

#define CTXTRUN2(ctxt1, ctxt2)                                                \
    if (dsp::ContextRun2<false, decltype(ctxt1), decltype(ctxt2)> contextRun{ \
            ctxt1, ctxt2})                                                    \
    contextRun = [&](auto ctxt1, auto ctxt2)
#define CTXTRUNREV2(ctxt1, ctxt2)                                            \
    if (dsp::ContextRun2<true, decltype(ctxt1), decltype(ctxt2)> contextRun{ \
            ctxt1, ctxt2})                                                   \
    contextRun = [&](auto ctxt1, auto ctxt2)

#define PROCESSBLOCK_                                \
    template <class Ctxt, class State>               \
    void processBlock(Ctxt ctxt, State &state) const \
    {                                                \
        CTXTRUN(ctxt) { process(ctxt, state); };     \
    }

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
