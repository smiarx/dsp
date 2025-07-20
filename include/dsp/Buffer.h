#pragma once

#include "Context.h"
#include "Utils.h"
#include <cstring>

namespace dsp
{
template <class T, std::size_t MinSize> class Buffer
{
    /* implements ring buffer */
  public:
    // offset to allow vectorized computing
    static constexpr auto kVecOffset = sizeof(batch<T>) / sizeof(T);
    static constexpr auto kBaseSize  = nextPow2(MinSize);
    static constexpr int kMask       = kBaseSize - 1;
    static constexpr auto kSize      = kBaseSize + kVecOffset;

    using Type = T;

    Buffer() = default;

    static constexpr auto getMinSize() { return MinSize; }

    void setData(T *data) { data_ = data; }
    T *getData() { return data_; }

    template <typename V, bool Vec = false, bool Safe = false>
    void write(int i, V val)
    {
        assert(data_ != nullptr);
        assert((0 <= i) && (i <= static_cast<int>(MinSize)));

        const auto pos = position(i);
        if constexpr (Vec) {
            storeBatch(data_[pos], val);
            /* copy end to beginning of buffer for vector continuity */
            if (static_cast<size_t>(pos) < V::kWidth) {
                // copy begining to end
                storeBatch(data_[kBaseSize], loadBatch(data_[0]));
            } else if (static_cast<size_t>(pos) > kBaseSize - V::kWidth) {
                // copy end to begining
                storeBatch(data_[0], loadBatch(data_[kBaseSize]));
            }
        } else {
            store(data_[pos], val);
            // safely write vor vec read even when writing scalars
            if constexpr (Safe && !std::is_same_v<T, batch<T>>) {
                if (static_cast<size_t>(pos) < kVecOffset) {
                    store(data_[kBaseSize + pos], val);
                }
            }
        }
    }
    template <typename V> void writeVec(int i, V val)
    {
        write<V, true>(i, val);
    }

    template <bool Vec = false> [[nodiscard]] auto read(int i) const
    {
        assert(data_ != nullptr);
        assert((0 <= i) && (i <= static_cast<int>(MinSize)));

        const auto pos = position(i);
        if constexpr (Vec) {
            return loadBatch(data_[pos]);
        } else {
            return load(data_[pos]);
        }
    }
    [[nodiscard]] auto readVec(int i) const { return read<true>(i); }

    void nextId(int incr) { id_ = (id_ + incr) & kMask; }

    // prepare for next block given ctxt
    template <class Ctxt> void nextBlock(Ctxt ctxt, bool checkLimits = false)
    {
        nextId(ctxt.getBlockSize());

        if (checkLimits && static_cast<size_t>(id_) < batch<T>::kWidth)
            // copy begining to end
            storeBatch(data_[kBaseSize], loadBatch(data_[0]));
    }

    [[nodiscard]] int getId() const { return id_; }

  private:
    int id_{kVecOffset};
    T *__restrict data_{nullptr};

    // give data position given i
    inline auto position(int i) const { return (id_ - i) & kMask; }
};

template <class T, size_t MinSize, bool Vec = false>
class BufferContext : public Context<T, Vec>
{
    using Parent  = Context<T, Vec>;
    using BufferT = Buffer<T, MinSize>;
    using SigType = typename Parent::SigType;

  public:
    BufferContext(T *in, int blockSize, const BufferT &buffer) :
        Parent(in, blockSize), buffer_(buffer)
    {
    }
    BufferContext(const Parent &ctxt, const BufferT &buffer) :
        Parent(ctxt), buffer_(buffer)
    {
    }

    [[nodiscard]] auto vec() const
    {
        return BufferContext<T, MinSize, true>(Context<T, Vec>::vec(), buffer_);
    }
    [[nodiscard]] auto scalar() const
    {
        return BufferContext<T, MinSize, false>(Context<T, Vec>::scalar(),
                                                buffer_);
    }

    BufferT &getBuffer() { return buffer_; }

    void next(int incr = Parent::kIncrSize)
    {
        buffer_.nextId(incr);
        Parent::next(incr);
    }

    void bufferSafeForVectorization() { buffer_.setLimits(); }

    [[nodiscard]] auto read(int i) const
    {
        return buffer_.template read<Vec>(i);
    }

    void write(int i, SigType val)
    {
        return buffer_.template write<SigType, Vec>(i, val);
    }

    void writeSafe(int i, SigType val)
    {
        return buffer_.template write<SigType, Vec, true>(i, val);
    }

  private:
    BufferT buffer_;
};
} // namespace dsp
