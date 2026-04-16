#pragma once

#include "dsp/analysis/ZeroCrossing.hpp"
#include "dsp/filters/IIRFilter.hpp"
#include "dsp/pitch/NSDF.hpp"
#include "dsp/pitch/PSOLA.hpp"
#include <dsp/core/Buffer.hpp>
#include <dsp/filters/AdaptiveFilter.hpp>
#include <dsp/filters/FIRFilter.hpp>
#include <dsp/filters/OnePole.hpp>
#include <dsp/pitch/PitchShift.hpp>

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

    void update(float forgetFactor, float warpIn, float warpOut, float ptop,
                float ptrho, float ptlambda, float shift);
    void prepare(float sampleRate);

  private:
    static constexpr auto kOrder = 8;

    dsp::RLS<float, kOrder> rls_{};
    dsp::WarpedAdaptiveFilter<float, kOrder> filter_{};
    decltype(filter_)::AnalyzeState warpAState_{};
    decltype(filter_)::ReconstructState warpRState_{};

    // amplitude follower
    dsp::OnePole<float> ampFollower_{};
    decltype(ampFollower_)::State ampFollowerState_{};
    decltype(ampFollower_)::State poleFollowerState_{};

    // zero crossing
    dsp::ZeroCrossing<float> zc_{};
    dsp::OnePole<float> zcOp_{};
    decltype(zcOp_)::State zcOpState_{};

    // pitch tracking
    bool chooseFastPitchTracker_{false};
    dsp::AdaptiveNotchFilter<float> fastPitchTracker_{0.999f, 0.99f};
    decltype(fastPitchTracker_)::State fastPitchTrackerState_{};
    dsp::NSDF<float, 300> pitchTracker_{};
    decltype(pitchTracker_)::DL<0> pitchTrackerDL_{};

    dsp::IIRFilter<float, 8> ptLowpass_{};
    decltype(ptLowpass_)::State ptLowpassState_{};
    dsp::OnePole<float> ptOp1_{};
    decltype(ptOp1_)::State ptOp1St_{};
    dsp::OnePole<float> ptOp2_{};
    decltype(ptOp2_)::State ptOp2St_{};

    // pitch shift
    float shift_{1};
    // dsp::PitchShift<float> pitchShift_{};
    // dsp::DelayLine<3000, nextTo(pitchTrackerDL_)> pitchShiftDL_{};
    dsp::PSOLA<float, 16> psola_{};
    dsp::DelayLine<8000, nextTo(pitchTrackerDL_)> pitchShiftDL_{};

    // resample
    static constexpr auto kDecFactor = 2;
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
