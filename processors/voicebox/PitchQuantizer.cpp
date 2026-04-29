#include "PitchQuantizer.h"
#include <dsp/core/Context.hpp>
#include <ranges>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

dsp::Pitch<float> PitchQuantizer::process(dsp::Pitch<float> pitch) noexcept
{
    auto octave        = pitch.octave();
    auto pitchInOctave = pitch - octave;

    auto quantizedPitches =
        quantizedPitches_ | std::views::take(nQuantizedPitches_);
    // chose pitch closest to our input pitch
    auto chosenPitch = *std::ranges::min_element(
        quantizedPitches, {}, [pitchInOctave](auto pitch) {
            return dsp::abs(pitch - pitchInOctave);
        });

    // put back in octave
    chosenPitch += octave;

    if (reset_) {
        reset_                = false;
        quantizedFilterState_ = chosenPitch.get();
        chosenPitch_          = chosenPitch;
    } else {
        // chosen pitch must have been chosen for some time
        if (chosenPitch != chosenPitch_) {
            if (chosenPitch != candidatePitch_) {
                candidatePitch_       = chosenPitch;
                candidatePitchLength_ = 0;
            }
            ++candidatePitchLength_;
            if (candidatePitchLength_ > confirmPitch_) {
                chosenPitch_ = chosenPitch;
            } else {
                chosenPitch = chosenPitch_;
            }
        }
        // smooth filter pitch
        auto chosenPitchValue = chosenPitch.get();
        quantizedFilter_.process(dsp::Context(&chosenPitchValue),
                                 quantizedFilterState_);
        chosenPitch = chosenPitchValue;
    }

    return chosenPitch;
};

void PitchQuantizer::setScale(dsp::Scale scale) noexcept
{
    nQuantizedPitches_ = 0;
    scale.forEachNote([this](dsp::Note note) {
        quantizedPitches_[nQuantizedPitches_] = dsp::Pitch<float>{note};
        ++nQuantizedPitches_;
    });
    // set first pitch an octave higher as last quantized value
    quantizedPitches_[nQuantizedPitches_] =
        quantizedPitches_[0] + dsp::Pitch<float>::octaveUnit();
    ++nQuantizedPitches_;
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
