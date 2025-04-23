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
    static constexpr auto kOffset   = In::kVectorSize;
    static constexpr auto kBaseSize = nextPow2(MinSize);
    static constexpr int kMask      = kBaseSize - 1;
    static constexpr auto kSize     = kBaseSize + kOffset;

    using Type = In;

    static constexpr auto getMinSize() { return MinSize; }

    void setBuffer(In *buffer) { buffer_ = buffer; }
    In *getBuffer() { return buffer_; }

    void write(int i, const In &x) { buffer_[(bufId_ - i) & kMask] = x; }
    template <bool Safe = true> void write(int i, const typename In::Vector &x)
    {
        const auto pos = (bufId_ - i) & kMask;
        if constexpr (Safe) {
            /* copy at the end of buffer for vector continuity */
            if (pos == 0) {
                buffer_[kBaseSize].toVector() = x;
            }
        }
        buffer_[pos].toVector() = x;
    }
    void write(int i, const typename In::Scalar &x)
    {
        buffer_[(bufId_ - i) & kMask].toScalar() = x;
    }
    [[nodiscard]] const In &read(int i) const
    {
        assert(i <= static_cast<int>(MinSize));
        return buffer_[(bufId_ - i) & kMask];
    }

    void nextBufId(int incr) { bufId_ = (bufId_ + incr) & kMask; }

    // prepare for next block given ctxt
    template <class Ctxt> void nextBlock(Ctxt ctxt)
    {
        nextBufId(ctxt.getBlockSize());
    }

    void setLimits() { buffer_[kBaseSize].toVector() = buffer_[0].toVector(); }

  private:
    bool isBufferOffset_{false};
    int bufId_{kOffset};
    In *__restrict buffer_{nullptr};
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

    [[nodiscard]] auto vec() const
    {
        return BufferContext<In, Buffer, true>(Context<In, Vectorize>::vec(),
                                               buffer_);
    }
    [[nodiscard]] auto scalar() const
    {
        return BufferContext<In, Buffer, false>(
            Context<In, Vectorize>::scalar(), buffer_);
    }

    Buffer &getBuffer() { return buffer_; }

    void next(int incr = Parent::kVecSize)
    {
        buffer_.nextBufId(incr);
        Parent::next(incr);
    }

    void bufferLimits() { buffer_.setLimits(); }

    [[nodiscard]] const auto &read(int i) const
    {
        auto &x = buffer_.read(i);
        return x.template toSignal<Vectorize>();
    }

  protected:
    Buffer buffer_;
};
} // namespace dsp
