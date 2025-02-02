#pragma once

#include <cstdint>

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

} // namespace dsp
