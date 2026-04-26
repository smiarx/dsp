#pragma once

#include <cassert>
#include <cstdint>
#include <dsp/cpu/defines.h>
#include <dsp/music/Note.hpp>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

class Scale
{
    // set a note representing a scale

    using bits = uint32_t;

  public:
    static constexpr auto kNumNotes = 12;
    static constexpr bits kMask     = (static_cast<bits>(1) << kNumNotes) - 1;

    constexpr explicit Scale(bits value) : value_(value) {}
    constexpr Scale(Note note) :
        value_(static_cast<bits>(1) << static_cast<bits>(note))
    {
    }

    constexpr bool operator==(const Scale &other) const = default;
    constexpr bool operator!=(const Scale &other) const = default;

    // combine scales
    constexpr Scale operator|(const Scale &other) const
    {
        return Scale{value_ | other.value_};
    }

    // combine scale and note
    constexpr Scale operator|(Note note) const { return *this | Scale(note); }

    constexpr Scale operator&(const Scale &other) const
    {
        return Scale{value_ & other.value_};
    }

    constexpr Scale operator^(const Scale &other) const
    {
        return Scale{value_ ^ other.value_};
    }

    constexpr Scale operator^(Note note) const { return *this ^ Scale(note); }

    constexpr Scale operator~() const { return Scale{(~value_) & kMask}; }

    // transpose scale: A becomes note
    constexpr Scale operator<<(Note note) const
    {
        auto noteNumber = static_cast<bits>(note);
        noteNumber %= kNumNotes;

        if (noteNumber == 0) return *this;

        auto result = (value_ << noteNumber) & kMask;
        result |= (value_ >> (kNumNotes - noteNumber)) & kMask;
        return Scale{result};
    }

    // transpose scale: note becomes A
    constexpr Scale operator>>(Note rhs) const
    {
        auto noteNumber = static_cast<bits>(rhs);
        noteNumber %= kNumNotes;

        if (noteNumber == 0) return *this;

        auto result = (value_ >> noteNumber) & kMask;
        result |= (value_ << (kNumNotes - noteNumber)) & kMask;
        return Scale{result};
    }

    // check if scale contains note
    [[nodiscard]] constexpr bool contains(Note note) const
    {
        return (value_ & Scale(note).value_) != 0;
    }

    // loop scale notes with function f
    template <typename F> constexpr void forEachNote(F f) const
    {
        for (std::size_t n = 0; n < kNumNotes; ++n) {
            auto note = static_cast<Note>(n);
            if (contains(note)) {
                f(note);
            }
        }
    }

  private:
    bits value_;
};

constexpr Scale operator|(Note lhs, Note rhs) { return Scale{lhs} | rhs; }
constexpr Scale operator~(Note note) { return ~Scale{note}; }

namespace scales
{
using enum Note;
static constexpr auto chromatic =
    A | As | B | C | Cs | D | Ds | E | F | Fs | G | Gs;

static constexpr auto minor           = A | B | C | D | E | F | G;
static constexpr auto minorHarmonic   = A | B | C | D | E | F | Gs;
static constexpr auto minorMelodic    = A | B | C | D | E | Fs | Gs;
static constexpr auto minorPentatonic = A | C | D | E | G;

static constexpr auto major           = (C | D | E | F | G | A | B) >> C;
static constexpr auto majorPentatonic = (C | D | E | G | A) >> C;

static constexpr auto ionian     = (A | B | C | D | E | F | G) >> C;
static constexpr auto dorian     = (A | B | C | D | E | F | G) >> D;
static constexpr auto phrygian   = (A | B | C | D | E | F | G) >> E;
static constexpr auto lydian     = (A | B | C | D | E | F | G) >> F;
static constexpr auto mixolydian = (A | B | C | D | E | F | G) >> G;
static constexpr auto aeolian    = (A | B | C | D | E | F | G) >> A;
static constexpr auto locrian    = (A | B | C | D | E | F | G) >> B;
} // namespace scales

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
