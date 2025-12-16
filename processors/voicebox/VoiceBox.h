#pragma once

#include "dsp/AdaptiveFilter.h"
#include "dsp/Buffer.h"
#include "dsp/FIRFilter.h"
#include "dsp/OnePole.h"
#include "dsp/PitchShift.h"

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

class VoiceBox
{
  public:
    VoiceBox()
    {
        buffer_.setData(bufferData_);
        bufferDec_.setData(bufferDecData_);
    }
    void process(const float *const *__restrict in,
                 float *const *__restrict out, int count);

    void update(float forgetFactor, float warpIn, float warpOut, float shift);
    void prepare(float sampleRate);

  private:
    static constexpr auto kOrder = 12;

    dsp::RLS<float, kOrder> rls_{};
    dsp::WarpedAdaptiveFilter<float, kOrder> filter_{};
    decltype(filter_)::AnalyzeState warpAState_{};
    decltype(filter_)::ReconstructState warpRState_{};

    // amplitude follower
    dsp::OnePole<float> ampFollower_{};
    decltype(ampFollower_)::State ampFollowerState_{};

    // pitch tracking
    dsp::AdaptiveNotchFilter<float> pitchTracker_{0.992f, 0.995f};
    decltype(pitchTracker_)::State pitchTrackerState_{};

    // pitch shift
    dsp::PitchShift<float> pitchShift_;
    dsp::DelayLine<3000> pitchShiftDL_;

    // resample
    static constexpr auto kDecFactor = 3;
    dsp::FIRDecimate<float, 21, kDecFactor> decimate_{};
    dsp::FIRInterpolate<float, 21, kDecFactor> interpolate_{};
    decltype(decimate_)::DL<0> decimateDL_{};
    decltype(interpolate_)::DL<nextTo(pitchShiftDL_)> interpolateDL_{};
    int decimateId_{};

    dsp::Buffer<float, nextTo(decimateDL_) + 512> buffer_;
    float bufferData_[decltype(buffer_)::kSize]{};
    dsp::Buffer<float, nextTo(interpolateDL_) + 512> bufferDec_;
    float bufferDecData_[decltype(bufferDec_)::kSize]{};
};
} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
