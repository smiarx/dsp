#pragma once

#include <cstdint>
#include <dsp/cpu/defines.h>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

enum class Note {
    A  = 0,
    As = 1,
    B  = 2,
    C  = 3,
    Cs = 4,
    D  = 5,
    Ds = 6,
    E  = 7,
    F  = 8,
    Fs = 9,
    G  = 10,
    Gs = 11,
};
static constexpr int kNumNotes = 12;

class Scale
{
  public:
    constexpr Scale(uint32_t value = (1 << 13) - 1) : value_(value) {}
    constexpr Scale(Note n) : value_(1 << static_cast<uint32_t>(n)) {}

    [[nodiscard]] constexpr auto value() const { return value_; }

    constexpr bool operator==(const Scale &other) const
    {
        return value_ == other.value_;
    }

    constexpr bool operator!=(const Scale &other) const
    {
        return value_ != other.value_;
    }

  private:
    uint32_t value_;
};

constexpr Scale operator|(Note lhs, Note rhs)
{
    return {(1u << static_cast<uint32_t>(lhs)) |
            (1u << static_cast<uint32_t>(rhs))};
}

constexpr Scale operator|(Scale lhs, Note rhs)
{
    return {lhs.value() | (1u << static_cast<uint32_t>(rhs))};
}

constexpr Scale operator|(Scale lhs, Scale rhs)
{
    return {lhs.value() | rhs.value()};
}

constexpr bool operator&(Scale lhs, Note rhs)
{
    return lhs.value() & (1u << static_cast<uint32_t>(rhs));
}

constexpr Scale operator<<(Scale lhs, Note rhs)
{
    auto noteNumber = static_cast<uint32_t>(rhs);
    auto value      = lhs.value();
    auto result     = (value << noteNumber) & ((1u << kNumNotes) - 1);
    result |= value >> (kNumNotes - noteNumber);
    return {result};
}

namespace scales
{
static constexpr auto minor =
    Note::A | Note::B | Note::C | Note::D | Note::E | Note::F | Note::G;
static constexpr auto major =
    Note::A | Note::B | Note::Cs | Note::D | Note::E | Note::Fs | Note::Gs;
} // namespace scales

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
