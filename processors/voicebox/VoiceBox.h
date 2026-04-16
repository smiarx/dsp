#pragma once

#include "PitchTracker.hpp"
#include "dsp/music/Pitch.hpp"
#include <atomic>
#include <dsp/analysis/ZFF.hpp>
#include <dsp/core/SampleRate.hpp>
#include <dsp/generators/GrainGenerator.hpp>
#include <dsp/music/Scale.hpp>
#include <dsp/pitch/PSOLA.hpp>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

using namespace dsp::units;

class VoiceBox
{
    using Hz      = dsp::freq::Hz<float>;
    using seconds = dsp::time::Seconds<float>;

  public:
    static constexpr auto kDefaultSampleRate  = dsp::SampleRate{48000._Hzf};
    static constexpr auto kMaxBlockSize       = 512;
    static constexpr auto kAnalysisWindowSize = 0.025_sf;
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

    void update(float shift, float formantShift, seconds correctSpeed,
                Hz a4Frequency);
    void prepare(Hz sampleRate);

    void setPitchShift(float shift);
    void setFormantShift(float shift);
    void setCorrectSpeed(seconds time);
    void setA4Frequency(Hz freq);
    void setScale(dsp::Scale scale);
    void setRange(PitchTracker::Range range);

    [[nodiscard]] float getPitchShift() const { return shift_; }
    [[nodiscard]] float getFormantShift() const { return formantShift_; }
    [[nodiscard]] float getCorrectSpeed() const { return correctSpeed_.get(); }
    [[nodiscard]] float getA4Frequency() const { return a4Frequency_.get(); }
    [[nodiscard]] auto getScale() const { return scale_; }
    [[nodiscard]] auto getRange() const { return pitchTracker_.getRange(); }

    [[nodiscard]] const std::atomic<float> &getAtomicFreq() const
    {
        return atomicFreq_;
    }

  private:
    float shift_{1.f};
    float formantShift_{1.f};
    seconds correctSpeed_{0.01_sf};
    Hz a4Frequency_{440._Hzf};
    dsp::Scale scale_{dsp::scales::chromatic};

    dsp::SampleRate<float> sampleRate_{kDefaultSampleRate};

    bool voiced_{};

    // preprocess
    dsp::IIRFilter<float, 8> dcblocker_{};
    decltype(dcblocker_)::State dcblockerState_{};

    // pitch tracking
    PitchTracker pitchTracker_{};

    dsp::OnePole<float> voicedSmooth_{};
    decltype(voicedSmooth_)::State voicedSmoothState_{};

    dsp::ZFF<float, kMaxAnalysisMarkers> zffAnalysis_{};
    dsp::DelayLine<1200> analysisDL_{};

    // compute pitch shift
    float computePitchShift(dsp::freq::NyquistRatio<float> inFreq);
    dsp::PitchReference<dsp::freq::NyquistRatio<float>> pitchRef_{};
    std::array<dsp::Pitch<float>, dsp::Scale::kNumNotes + 1> scalePitches_{};
    size_t nScalePitches_{};
    dsp::OnePole<float> correctionFilter_{};
    decltype(correctionFilter_)::State correctionFilterState_{};

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
