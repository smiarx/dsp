#pragma once

#include <cmath>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T> static constexpr T nextAlignedOffset(T off, T align)
{
    off += align - 1;
    return off - (off & (align - 1));
}

static constexpr std::size_t nextPow2(std::size_t size)
{
    size -= 1;
    size |= size >> 1;
    size |= size >> 2;
    size |= size >> 4;
    size |= size >> 8;
    size |= size >> 16;
    if constexpr (sizeof(std::size_t) > 4) {
        size |= size >> 32;
    }
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
    return pow(F(10), db / F(20));
}

template <typename T, typename U>
static constexpr auto expScale(T min, T max, U x)
{
    return min * std::pow(max / min, x);
}

////////////////////////////////////////////////////////////////////
/// Modular arithmetic
template <int Mod> static bool modCmpGE(int id1, int id2)
{
    static_assert(isPow2(Mod));
    // check if id1 >= id2 modulo [kBaseSize]
    // (e.g 5 [2048] >= 2030 [2048])
    constexpr auto kMidSize = Mod >> 1;
    return ((id1 - id2) & kMidSize) == 0;
}
template <int Mod> static bool modCmpGT(int id1, int id2)
{
    static_assert(isPow2(Mod));
    // check if id1 > id2 modulo [kBaseSize]
    // (e.g 5 [2048] > 2030 [2048])
    return modCmpGE<Mod>(id1, id2 + 1);
}
template <int Mod> static bool modCmpEq(int id1, int id2)
{
    static_assert(isPow2(Mod));
    // check if id1 == id2 modulo [kBaseSize]
    // (e.g 5 [2048] > 2053 [2048])
    return ((id1 - id2) & (Mod - 1)) == 0;
}

template <int Mod, typename Id1, typename Id2>
static auto modAdd(Id1 id1, Id2 id2)
{
    static_assert(isPow2(Mod));
    if constexpr (std::is_floating_point_v<Id2>) {
        constexpr auto kMod     = static_cast<Id2>(Mod);
        constexpr auto kMidSize = static_cast<Id2>(Mod >> 1);
        auto add                = static_cast<Id2>(id1) + id2;
        if (add > kMidSize) {
            return -kMod + add;
        } else if (add < -kMidSize) {
            return kMod + add;
        }
        return add;
    } else {
        return (id1 + id2) & (Mod - 1);
    }
}

template <int Mod, typename Id1, typename Id2>
static auto modDiff(Id1 id1, Id2 id2)
{
    static_assert(isPow2(Mod));
    if constexpr (std::is_floating_point_v<Id2>) {
        return modAdd<Mod>(id1, -id2);
    } else {
        constexpr auto kMidSize = Mod >> 1;
        constexpr auto kMidMask = kMidSize - 1;
        auto diff               = (id1 - id2);
        return (diff & kMidMask) - (diff & kMidSize);
    }
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
