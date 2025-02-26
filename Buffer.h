#pragma once

#include "Context.h"
#include "Utils.h"
#include <cstring>

namespace dsp
{
template <class In, std::size_t MinSize = 0> class Buffer
{
  public:
    /* let the Offset first element be free and copy and of buffer to it so that
     * we can always retrieve Vec */
    static constexpr auto Offset   = In::VectorSize;
    static constexpr auto BaseSize = nextPow2(MinSize);
    static constexpr auto Mask     = BaseSize - 1;
    static constexpr auto Size     = BaseSize + Offset;

    using Type = In;

    static constexpr auto getMinSize() { return MinSize; }

    void setBuffer(In *buffer) { buffer_ = buffer; }
    In *getBuffer() { return buffer_; }

    void write(int i, const In &x) { buffer_[(bufId_ - i) & Mask] = x; }
    template <bool Safe = true> void write(int i, const typename In::Vector &x)
    {
        const auto pos = (bufId_ - i) & Mask;
        if constexpr (Safe) {
            /* copy at the end of buffer for vector continuity */
            if (pos == 0) {
                buffer_[BaseSize].toVector() = x;
            }
        }
        buffer_[pos].toVector() = x;
    }
    void write(int i, const typename In::Scalar &x)
    {
        buffer_[(bufId_ - i) & Mask].toScalar() = x;
    }
    const In &read(size_t i) const
    {
        assert(i <= MinSize);
        return buffer_[(bufId_ - i) & Mask];
    }

    void nextBufId(int incr) { bufId_ = (bufId_ + incr) & Mask; }

    void setLimits() { buffer_[BaseSize].toVector() = buffer_[0].toVector(); }

  private:
    bool isBufferOffset_{false};
    int bufId_{Offset};
    In *buffer_;
};

template <class In, class Buffer, bool Vectorize = false>
class BufferContext : public Context<In, Vectorize>
{
    using Parent = Context<In, Vectorize>;

  public:
    BufferContext(In *in, int blockSize, const Buffer &buffer) :
        Parent(in, blockSize), buffer_(buffer)
    {
    }
    BufferContext(const Parent &ctxt, const Buffer &buffer) :
        Parent(ctxt), buffer_(buffer)
    {
    }

    auto vec() const
    {
        return BufferContext<In, Buffer, true>(Context<In, Vectorize>::vec(),
                                               buffer_);
    }
    auto scalar() const
    {
        return BufferContext<In, Buffer, false>(
            Context<In, Vectorize>::scalar(), buffer_);
    }

    Buffer &getBuffer() { return buffer_; }

    void next(int incr = Parent::VecSize)
    {
        buffer_.nextBufId(incr);
        Parent::next(incr);
    }

    void nextBlock()
    {
        const auto blockSize = Parent::getBlockSize();
        next(blockSize);
    }

    void bufferLimits() { buffer_.setLimits(); }

    const auto &read(int i) const
    {
        auto &x = buffer_.read(i);
        return x.template toSignal<Vectorize>();
    }

    void save(Buffer &buffer) { buffer = buffer_; }

  protected:
    Buffer buffer_;
};
} // namespace dsp
