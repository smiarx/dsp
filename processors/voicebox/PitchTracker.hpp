#pragma once

#include "dsp/filters/Amplitude.hpp"
#include <dsp/analysis/ZeroCrossing.hpp>
#include <dsp/core/Buffer.hpp>
#include <dsp/filters/AdaptiveFilter.hpp>
#include <dsp/filters/IIRFilter.hpp>
#include <dsp/filters/OnePole.hpp>
#include <dsp/filters/VAFilters.hpp>
#include <dsp/math/HMM.hpp>
#include <dsp/pitch/NSDF.hpp>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

class PitchTracker
{
    static constexpr auto kMaxBlockSize = 512;
    static constexpr auto kDecFactor    = 8;

    static constexpr auto kDefaultSampleRate    = 48000.f;
    static constexpr auto kTargetDownSampleRate = 6000.f;
    static constexpr auto kMaxDetectedFreq      = 1200.f;

    static constexpr auto kVoicedProbabilityStartTracking  = 0.8f;
    static constexpr auto kVoicedProbabilityStopConfirm    = 0.2f;
    static constexpr auto kAutoCorrelationPeakStopTracking = 0.f;

  public:
    PitchTracker()
    {
        buffer_.setData(bufferData_.data());
        bufferDec_.setData(bufferDecData_.data());
    };

    void process(const float *const *__restrict in,
                 float *const *__restrict out, int count);

    void prepare(float sampleRate);

    enum class Range {
        kLow,
        kMedium,
        kHigh,
    };

    void setRange(Range range);
    [[nodiscard]] Range getRange() const noexcept { return range_; }

  private:
    Range range_{Range::kHigh};

    dsp::interp::PolyPoint<float> findNSDFPeriod();
    float smoothVoicedProbability(float p);

    dsp::va::OnePole<float, dsp::va::kLowPass> softLowpass_{};
    decltype(softLowpass_)::State softLowpassState_{};

    // decimate filter
    dsp::IIRFilter<float, 8> nyquistFilter_{};
    decltype(nyquistFilter_)::State nyquistFilterState_{};

    // decimate
    int decimateId_{};

    float downPeriod_{};
    float downProb_{1.f};
    float period_{};

    // nsdf
    dsp::NSDF<float, 120> nsdf_{};
    decltype(nsdf_)::DL<0> nsdfDL_{};
    float autoCorrelationThreshold_{0.84f};

    enum HarmonicJump {
        kJumpOctaveDown,
        kNoJump,
        kJumpOctaveUp,
        kJump3rdHarmonicUp,
        kNumJumps,
    } harmonicJumpType_{kNoJump};
    int harmonicJumpCounter_{};
    void resetHarmonicJump();

    // voiced detection
    dsp::Amplitude<float> ampFollower_{};
    decltype(ampFollower_)::State ampState_{};
    decltype(ampFollower_)::State ampFilteredState_{};
    // zero crossing voicing detection
    dsp::ZeroCrossing<float> zeroCrossing_{};
    dsp::Amplitude<float> zeroCrossingFollower_{};
    decltype(zeroCrossingFollower_)::State zeroCrossingFollowerState_{};

    dsp::HMM2<float> hmm2_{0.92f};

    // track
    dsp::NSDFTrack<float, 16, 960> nsdfTrack_{};
    decltype(nsdfTrack_)::DL<0> nsdfTrackDL_{};

    dsp::va::SVF<float, dsp::va::kLowPass> smoothFreq_{};
    decltype(smoothFreq_)::State smoothFreqState_{};

    bool voiced_{};
    int nConfirmPitch_{};
    float confirmPeriodMult_{3.f};

    float freqScale_;
    float decFactor_{kDefaultSampleRate / kTargetDownSampleRate};
    int iDecFactor_{static_cast<int>(decFactor_)};
    float minPeriodDetected_{kTargetDownSampleRate / kMaxDetectedFreq};

    dsp::Buffer<float, nextTo(nsdfDL_) + kMaxBlockSize> bufferDec_{};
    std::array<float, decltype(bufferDec_)::kSize> bufferDecData_{};

    dsp::Buffer<float, nextTo(nsdfTrackDL_) + kMaxBlockSize> buffer_;
    std::array<float, decltype(buffer_)::kSize> bufferData_{};
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
