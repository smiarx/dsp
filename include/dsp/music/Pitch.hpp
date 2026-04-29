#pragma once

#include <dsp/core/Units.hpp>
#include <dsp/music/Scale.hpp>

namespace dsp
{

inline namespace DSP_ARCH_NAMESPACE
{

// pitch octave based
template <typename T> struct Pitch : public Unit<Pitch, T> {

    static consteval Pitch octaveUnit() { return {T(1)}; }

    constexpr Pitch(T value = T(0)) : Unit<Pitch, T>{value} {}

    explicit constexpr Pitch(Note note, Note baseNote = Note::A) :
        Unit<Pitch, T>(noteToPitch(note, baseNote))
    {
    }

    // Converts pitch in octaves to a frequency ratio.
    [[nodiscard]] T toRatio() const { return std::exp2(this->get()); }

    [[nodiscard]] Pitch octave() const { return {dsp::floor(this->get())}; }

  private:
    constexpr T noteToPitch(Note note, Note baseNote = Note::A)
    {
        auto n = static_cast<int>(note) - static_cast<int>(baseNote);
        n      = (n + Scale::kNumNotes) % Scale::kNumNotes;
        return static_cast<T>(n) / static_cast<T>(Scale::kNumNotes);
    }
};

template <freq::Frequency F> class PitchReference
{
    using value_type = typename F::type;

  public:
    PitchReference() = default;
    explicit PitchReference(F freq) : reference_(std::log2(freq.get())) {}

    // return freq as pitch where pitch == 0 if freq == reference, every step of
    // 1 represents an octave
    auto toPitch(F freq) { return Pitch{std::log2(freq.get())} - reference_; }

    F toFreq(Pitch<value_type> pitch)
    {
        return {(pitch + reference_).toRatio()};
    }

  private:
    Pitch<value_type> reference_{};
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
