#pragma once

#include "Context.h"
#include <array>

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
    static constexpr auto Size   = nextPow2(MinSize);
    static constexpr auto Mask   = Size - 1;

    using Type = In;

    void setBuffer(In *buffer) { memset(buffer, 0, MinSize*sizeof(In)); buffer_ = buffer;}
    In *getBuffer() { return buffer_; }


    void write(int i, const In &x) { buffer_[(bufId_ + i) & Mask] = x; }
    const In &read(int i) const { return buffer_[(bufId_ + i) & Mask]; }

    template<int S> std::array<In,S> readContigous(int i) const {
        std::array<In,S> val;
        auto* ptr = &buffer_[bufId_+i];
        if(bufId_ + i + S < Size)
        {
            std::copy(ptr,ptr+S,val.begin());
        }
        else
        {
            auto dif = Size-bufId_;
            std::copy(ptr,ptr+dif,val.begin());
            ptr = buffer_;
            std::copy(ptr,ptr+(S-dif),val.begin()+dif);
        }
        return val;
    }

    void nextBufId(int incr) { bufId_ = (bufId_ - incr) & Mask; }

  private:
    int bufId_{0};
    In *buffer_;
};

template <class In, class Buffer> class BufferContext : public Context<In>
{
  public:
    BufferContext(In *in, Buffer &buffer) : Context<In>(in), buffer_(buffer) {}
    BufferContext(BufferContext &ctxt, int step = 1) :
        Context<In>(ctxt, step), buffer_(ctxt.buffer_)
    {
    }

    Buffer &getBuffer() { return buffer_; }
    void nextBufId(int incr) { buffer_.nextBufId(incr); }

    void next(int incr = 1)
    {
        nextBufId(incr);
        Context<In>::next(incr*Context<In>::getStep());
    }
    void nextBlock() { next(Context<In>::getBlockSize()); }

    void save(Buffer &buffer) { buffer = buffer_; }

  protected:
    Buffer buffer_;
};
} // namespace dsp
