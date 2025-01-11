#pragma once

namespace dsp
{
static constexpr int nextAlignedOffset(int off, int align)
{
    off += align - 1;
    return off - off % align;
}

} // namespace dsp
