#pragma once

#include "PitchTracker.hpp"
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
    static constexpr auto kAnalysisWindowSize = 0.015f;
    static constexpr auto kMaxAnalysisMarkers = 12;
    VoiceBox() { buffer_.setData(bufferData_); }
    void process(const float *const *__restrict in,
                 float *const *__restrict out, int count);

    void update(float shift);
    void prepare(float sampleRate);

  private:
    bool voiced_{};
    float sampleRate_{};

    // preprocess
    dsp::va::OnePole<float, dsp::va::kHighPass> dcblocker_{};
    decltype(dcblocker_)::State dcblockerState_{};

    // pitch tracking
    PitchTracker pitchTracker_{};

    dsp::ZFF<float, kMaxAnalysisMarkers> zffAnalysis_{};
    dsp::DelayLine<1200> analysisDL_{};

    float shift_{1.f};
    dsp::PSOLA<float, 16> psola_{};
    dsp::DelayLine<10000, nextTo(analysisDL_)> synthDL_{};

    dsp::Buffer<float, nextTo(synthDL_) + kMaxBlockSize> buffer_{};
    float bufferData_[decltype(buffer_)::kSize]{};

    float pitchBuffer_[kMaxBlockSize]{};
};
} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
