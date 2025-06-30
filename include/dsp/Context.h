#pragma once

#include <cassert>
#include <cstdlib>
#include <tuple>

template <typename T> T load(T x) { return x; }
template <typename T> void store(T &dest, T val) { dest = val; }

namespace dsp
{
template <typename T> class Context
{
  public:
    Context(T *data, int blockSize = 1) : blockSize_(blockSize), data_(data) {}
    Context(const Context &ctxt)            = default;
    Context &operator=(const Context &ctxt) = default;

    using Type = T;

    auto getInput() const { return load(*data_); }
    template <typename S> void setOutput(S value) { store(*data_, value); }

    void next(int incr = 1) { nextData(incr); }

    [[nodiscard]] int getBlockSize() const { return blockSize_; }
    void setBlockSize(int blockSize) { blockSize_ = blockSize; }

    void setData(T *data) { data_ = data; }
    [[nodiscard]] const T *getData() const { return data_; }

  protected:
    void nextData(int incr) { data_ += incr; }

  private:
    int blockSize_;
    T *__restrict data_;
};

#define contextFor(ctxt)                                          \
    for (auto [c, n] = std::tuple{ctxt, 0}; n < c.getBlockSize(); \
         n += 1, c.next())

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
