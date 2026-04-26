#pragma once

#include "PitchTracker.hpp"
#include "dsp/generators/Noise.hpp"
#include <atomic>
#include <dsp/analysis/ZFF.hpp>
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

    VoiceBox()
    {
        buffer_.setData(bufferData_);
        rls_.setForgetFactor(0.9995f);
    }

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

    dsp::ZFF<float, kMaxAnalysisMarkers> zffAnalysis_{};
    dsp::DelayLine<1200> analysisDL_{};

    // RLS
    static constexpr auto kOrder = 8;
    dsp::Noise<float> noise_{};
    dsp::RLS<float, kOrder> rls_{};
    dsp::WarpedAdaptiveFilter<float, kOrder> filter_{};
    decltype(filter_)::AnalyzeState warpAState1_{};
    decltype(filter_)::AnalyzeState warpAState_{};
    decltype(filter_)::ReconstructState warpRState_{};

    float shift_{1.f};
    float formantShift_{1.f};
    dsp::PSOLA<float, 16> psola_{};
    dsp::DelayLine<10000, nextTo(analysisDL_)> synthDL_{};

    dsp::Buffer<float, nextTo(synthDL_) + kMaxBlockSize> buffer_{};
    float bufferData_[decltype(buffer_)::kSize]{};

    float pitchBuffer_[kMaxBlockSize]{};
    float voicedBuffer_[kMaxBlockSize]{};
    std::atomic<float> atomicFreq_;
};
} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
