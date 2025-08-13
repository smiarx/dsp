#pragma once

#include "dsp/AdaptiveFilter.h"
#include "dsp/Buffer.h"
#include "dsp/FIRFilter.h"

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

    void update(float forgetFactor, float warpIn, float warpOut);

  private:
    static constexpr auto kOrder = 12;

    dsp::RLS<float, kOrder> rls_{};
    dsp::WarpedAdaptiveFilter<float, kOrder> filter_{};
    decltype(filter_)::AnalyzeState warpAState_{};
    decltype(filter_)::ReconstructState warpRState_{};

    // resample
    dsp::FIRDecimate<float, 21, 3> decimate_{};
    dsp::FIRInterpolate<float, 21, 3> interpolate_{};
    decltype(decimate_)::DL<0> decimateDL_{};
    decltype(interpolate_)::DL<0> interpolateDL_{};
    int decimateId_{};

    dsp::Buffer<float, nextTo(decimateDL_) + 512> buffer_;
    float bufferData_[decltype(buffer_)::kSize]{};
    dsp::Buffer<float, nextTo(interpolateDL_) + 512> bufferDec_;
    float bufferDecData_[decltype(bufferDec_)::kSize]{};
};
} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
