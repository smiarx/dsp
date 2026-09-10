#pragma once

#include <dsp/core/Units.hpp>
#include <dsp/filters/OnePole.hpp>
#include <dsp/filters/VAFilters.hpp>
#include <dsp/music/Pitch.hpp>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

class PitchQuantizer
{

  public:
    void setScale(dsp::Scale scale) noexcept;

    void setQuantizeSpeed(dsp::time::Samples<float> samples) noexcept
    {
        constexpr auto kConvRatio = 0.05f;
        quantizedFilter_.setConvergenceTime(kConvRatio, samples);
    }

    void setQuantizeConfirmTime(dsp::time::Samples<float> samples) noexcept
    {
        confirmPitch_ = static_cast<size_t>(samples.get());
    }

    void setVariationHpCutoff(dsp::freq::NyquistRatio<float> freq) noexcept
    {
        hpfilter_.setFreq(freq);
    }

    void setVariationAmount(float variationAmount) noexcept
    {
        variationAmount_ = variationAmount;
    }
    [[nodiscard]] float getVariationAmount() const noexcept
    {
        return variationAmount_;
    }

    dsp::Pitch<float> process(dsp::Pitch<float> pitch) noexcept;

    void reset() noexcept { reset_ = true; }

  private:
    // pitch candidate
    dsp::Pitch<float> chosenPitch_{};
    dsp::Pitch<float> candidatePitch_{};
    size_t candidatePitchLength_{};
    size_t confirmPitch_{};

    // array of quantized pitches
    std::array<dsp::Pitch<float>, dsp::Scale::kNumNotes + 1>
        quantizedPitches_{};
    size_t nQuantizedPitches_{};

    // filter quantized pitch
    dsp::OnePole<float> quantizedFilter_{};
    decltype(quantizedFilter_)::State quantizedFilterState_{};

    // variations filter
    dsp::va::OnePole<float, dsp::va::kHighPass> hpfilter_{};
    decltype(hpfilter_)::State hpfilterState_{};
    float variationAmount_{};

    // reset quantization
    bool reset_{};
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
