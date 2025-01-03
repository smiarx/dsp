#pragma once

#include "Context.h"

namespace dsp
{
template <class In, std::size_t MinSize = 0> class Buffer
{
    static constexpr std::size_t nextPow2(std::size_t size)
    {
        size -= 1;
        size |= size >> 1;
        size |= size >> 2;
        size |= size >> 4;
        size |= size >> 8;
        size |= size >> 16;
        size |= size >> 32;
        size += 1;
        return size;
    }

  public:
    static constexpr auto Offset = MinSize;
    static constexpr auto Size = nextPow2(MinSize);
    static constexpr auto Mask = Size - 1;

    using Type = In;

    void setBuffer(In *buffer) { buffer_ = buffer; }
    const In &read(int i) const { return buffer_[(bufId_ + i) & Mask]; }
    void write(int i, const In &x) { buffer_[(bufId_ + i) & Mask] = x; }

    void nextBufId(int incr) { bufId_ = (bufId_ - incr) & Mask; }

  private:
    int bufId_{0};
    In *buffer_;
};

template <class In, class Buffer> class BufferContext : public Context<In>
{
  public:
    BufferContext(In *in, Buffer &buffer) : Context<In>(in), buffer_(buffer) {}

    Buffer &getBuffer() { return buffer_; }
    void nextBufId(int incr) { buffer_.nextBufId; }

    void next(int incr = 1)
    {
        buffer_.nextBufId(incr);
        Context<In>::next(incr);
    }

  protected:
    Buffer buffer_;
};
} // namespace dsp
