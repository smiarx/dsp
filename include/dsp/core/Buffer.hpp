#pragma once

#include <cstring>
#include <dsp/core/Context.hpp>
#include <dsp/utils/Utils.hpp>

namespace dsp
{
using bufferId_t = int;

inline namespace DSP_ARCH_NAMESPACE
{
template <class T, std::size_t MinSize> class Buffer
{
    /* implements ring buffer */
  public:
    // offset to allow vectorized computing
    static constexpr auto kVecOffset  = sizeof(batch<T>) / sizeof(T);
    static constexpr auto kBaseSize   = nextPow2(MinSize);
    static constexpr bufferId_t kMask = kBaseSize - 1;
    static constexpr auto kSize       = kBaseSize + kVecOffset;

    using Type = T;

    Buffer() = default;

    static constexpr auto getMinSize() { return MinSize; }

    void setData(T *data) { data_ = data; }
    T *getData() { return data_; }

    template <typename V, bool Vec = false, bool Safe = false>
    void write(bufferId_t i, V val)
    {
        assert(data_ != nullptr);
        assert((0 <= i) && (i <= static_cast<bufferId_t>(MinSize)));

        const auto pos = position(i);
        if constexpr (Vec) {
            storeBatch(data_[pos], val);
            /* copy end to beginning of buffer for vector continuity */
            if (static_cast<size_t>(pos) < V::kWidth) {
                // copy begining to end
                dsp::DSP_ARCH_NAMESPACE::storeBatch(
                    data_[kBaseSize],
                    dsp::DSP_ARCH_NAMESPACE::loadBatch(data_[0]));
            } else if (static_cast<size_t>(pos) > kBaseSize - V::kWidth) {
                // copy end to begining
                dsp::DSP_ARCH_NAMESPACE::storeBatch(
                    data_[0],
                    dsp::DSP_ARCH_NAMESPACE::loadBatch(data_[kBaseSize]));
            }
        } else {
            store(data_[pos], val);
            // safely write vor vec read even when writing scalars
            if constexpr (Safe && !std::is_same_v<T, batch<T>>) {
                if (static_cast<size_t>(pos) < kVecOffset) {
                    dsp::DSP_ARCH_NAMESPACE::store(
                        data_[static_cast<int>(kBaseSize) + pos], val);
                }
            }
        }
    }
    template <typename V> void writeVec(bufferId_t i, V val)
    {
        write<V, true>(i, val);
    }

    template <bool Vec = false> [[nodiscard]] auto read(bufferId_t i) const
    {
        assert(data_ != nullptr);
        assert((0 <= i) && (i <= static_cast<bufferId_t>(MinSize)));

        const auto pos = position(i);
        if constexpr (Vec) {
            return dsp::DSP_ARCH_NAMESPACE::loadBatch(data_[pos]);
        } else {
            return dsp::DSP_ARCH_NAMESPACE::load(data_[pos]);
        }
    }
    [[nodiscard]] auto readVec(bufferId_t i) const { return read<true>(i); }

    void nextId(bufferId_t incr) { id_ = (id_ + incr) & kMask; }

    // prepare for next block given ctxt
    template <class Ctxt> void nextBlock(Ctxt ctxt, bool checkLimits = false)
    {
        nextId(ctxt.getBlockSize());

        if (checkLimits && static_cast<size_t>(id_) < batch<T>::kWidth)
            // copy begining to end
            dsp::DSP_ARCH_NAMESPACE::storeBatch(
                data_[kBaseSize], dsp::DSP_ARCH_NAMESPACE::loadBatch(data_[0]));
    }

    [[nodiscard]] bufferId_t getId() const { return id_; }

  private:
    bufferId_t id_{kVecOffset};
    T *__restrict data_{nullptr};

    // give data position given i
    inline auto position(bufferId_t i) const { return (id_ - i) & kMask; }
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

    void next(bufferId_t incr = Parent::kIncrSize)
    {
        buffer_.nextId(incr);
        Parent::next(incr);
    }

    void next(int incr, bufferId_t bufferIncr)
    {
        buffer_.nextId(bufferIncr);
        Parent::next(incr);
    }

    void bufferSafeForVectorization() { buffer_.setLimits(); }

    [[nodiscard]] auto read(bufferId_t i) const
    {
        return buffer_.template read<Vec>(i);
    }

    void write(bufferId_t i, SigType val)
    {
        return buffer_.template write<SigType, Vec>(i, val);
    }

    void writeSafe(bufferId_t i, SigType val)
    {
        return buffer_.template write<SigType, Vec, true>(i, val);
    }

    [[nodiscard]] auto getId() const { return buffer_.getId(); }

    // id greater than
    bool idGreaterThan(bufferId_t otherId)
    {
        return modCmpGT<BufferT::kBaseSize>(getId(), otherId);
    }
    // id small than
    bool idLessThan(bufferId_t otherId)
    {
        return modCmpGT<BufferT::kBaseSize>(otherId, getId());
    }

    // id greater than
    bool idGreaterEq(bufferId_t otherId)
    {
        return modCmpGE<BufferT::kBaseSize>(getId(), otherId);
    }
    // id small than
    bool idLessEq(bufferId_t otherId)
    {
        return modCmpGE<BufferT::kBaseSize>(otherId, getId());
    }

    bool idEqual(bufferId_t otherId)
    {
        return modCmpEq<BufferT::kBaseSize>(getId(), otherId);
    }

    template <typename IdType> auto idDiff(IdType otherId)
    {
        return dsp::modDiff<BufferT::kBaseSize>(getId(), otherId);
    }

    template <typename IdType> auto modDiff(IdType id1, IdType id2)
    {
        return dsp::modDiff<BufferT::kBaseSize>(id1, id2);
    }

    template <typename IdType> auto modAdd(IdType id1, IdType id2)
    {
        return dsp::modAdd<BufferT::kBaseSize>(id1, id2);
    }

  private:
    BufferT buffer_;
};
} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
