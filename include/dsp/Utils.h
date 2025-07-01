#pragma once

#include <cmath>

namespace dsp
{
static constexpr int nextAlignedOffset(int off, int align)
{
    off += align - 1;
    return off - off % align;
}

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

static constexpr bool isPow2(int v) { return v && ((v & (v - 1)) == 0); }

static constexpr unsigned int ilog2(unsigned int v)
{
    unsigned int c = 32; // c will be the number of zero bits on the right
    v &= static_cast<unsigned int>(-signed(v));
    if (v) c--;
    if (v & 0x0000FFFF) c -= 16;
    if (v & 0x00FF00FF) c -= 8;
    if (v & 0x0F0F0F0F) c -= 4;
    if (v & 0x33333333) c -= 2;
    if (v & 0x55555555) c -= 1;
    return c;
}

template <typename F> static constexpr auto db2gain(F db)
{
    return std::pow(F(10), db / F(20));
}

template <typename F> static constexpr auto expScale(F min, F max, F x)
{
    return min * std::pow(max / min, x);
}

} // namespace dsp
