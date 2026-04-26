#pragma once

#include "PitchTracker.hpp"
#include <atomic>
#include <dsp/analysis/ZFF.hpp>
#include <dsp/generators/GrainGenerator.hpp>
#include <dsp/pitch/PSOLA.hpp>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

class VoiceBox
{
  public:
    static constexpr auto kMaxBlockSize       = 512;
    static constexpr auto kAnalysisWindowSize = 0.025f;
    static constexpr auto kMaxAnalysisMarkers = kMaxBlockSize;

    static constexpr auto kMarkerVoicedProbUpdate       = 0.9f;
    static constexpr auto kMarkerDetectedPitchRatio     = 0.8f;
    static constexpr auto kMarkerSimilarPitchLowFactor  = 0.9f;
    static constexpr auto kMarkerSimilarPitchHighFactor = 1.1f;

    static constexpr auto kPsolaHysteresisStart    = 0.2f;
    static constexpr auto kPsolaHysteresisStop     = 0.1f;
    static constexpr auto kGrainGenHysteresisStart = 0.8f;
    static constexpr auto kGrainGenHysteresisStop  = 0.95f;

    VoiceBox() { buffer_.setData(bufferData_); }

    void process(const float *const *__restrict in,
                 float *const *__restrict out, int count);

    void update(float shift, float formantShift);
    void prepare(float sampleRate);

    void setPitchShift(float shift);
    void setFormantShift(float shift);
    void setRange(PitchTracker::Range range);

    [[nodiscard]] float getPitchShift() const { return shift_; }
    [[nodiscard]] float getFormantShift() const { return formantShift_; }
    [[nodiscard]] auto getRange() const { return pitchTracker_.getRange(); }

    [[nodiscard]] const std::atomic<float> &getAtomicFreq() const
    {
        return atomicFreq_;
    }

  private:
    bool voiced_{};
    float sampleRate_{};

    // preprocess
    dsp::IIRFilter<float, 8> dcblocker_{};
    decltype(dcblocker_)::State dcblockerState_{};

    // pitch tracking
    PitchTracker pitchTracker_{};

    dsp::OnePole<float> voicedSmooth_{};
    decltype(voicedSmooth_)::State voicedSmoothState_{};

    dsp::ZFF<float, kMaxAnalysisMarkers> zffAnalysis_{};
    dsp::DelayLine<1200> analysisDL_{};

    float shift_{1.f};
    float formantShift_{1.f};
    dsp::PSOLA<float, 16> psola_{};
    dsp::DelayLine<10000, nextTo(analysisDL_)> synthDL_{};
    dsp::GrainGenerator<float, 16> grainGenerator_{};

    dsp::Buffer<float, nextTo(synthDL_) + kMaxBlockSize> buffer_{};
    float bufferData_[decltype(buffer_)::kSize]{};

    float pitchBuffer_[kMaxBlockSize]{};
    float voicedBuffer_[kMaxBlockSize]{};
    std::atomic<float> atomicFreq_;
};
} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
