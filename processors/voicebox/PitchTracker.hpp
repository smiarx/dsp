#pragma once

#include "dsp/analysis/ZeroCrossing.hpp"
#include "dsp/filters/Notch.hpp"
#include <dsp/core/Buffer.hpp>
#include <dsp/filters/AdaptiveFilter.hpp>
#include <dsp/filters/IIRFilter.hpp>
#include <dsp/filters/OnePole.hpp>
#include <dsp/filters/VAFilters.hpp>
#include <dsp/pitch/NSDF.hpp>
#include <dsp/math/HMM.hpp>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

class PitchTracker
{
    static constexpr auto kMaxBlockSize       = 512;
    static constexpr auto kDecFactor = 8;

  public:
    PitchTracker() { bufferDec_.setData(bufferDecData_.data()); };

    void process(const float *const *__restrict in,
                 float *const *__restrict out, int count);

    void prepare(float sampleRate);

    float findNSDFPeriod(float zeroCrossing);

  private:
    dsp::va::OnePole<float, dsp::va::kLowPass> softLowpass_{};
    decltype(softLowpass_)::State softLowpassState_{};

    // decimate filter
    dsp::IIRFilter<float, 8> nyquistFilter_{};
    decltype(nyquistFilter_)::State nyquistFilterState_{};

    // decimate
    int decimateId_{};

    // nsdf
    dsp::NSDF<float, 120> nsdf_{};
    decltype(nsdf_)::DL<0> nsdfDL_{};

    // zero crossing voicing detection
    dsp::ZeroCrossing<float> zeroCrossing_{};
    dsp::OnePole<float> zeroCrossingFollower_{};
    decltype(zeroCrossingFollower_)::State zeroCrossingFollowerState_{};
    float zeroCrossingBuffer_[kMaxBlockSize]{};
    // voicing detection
    dsp::HMM2<float> hmm2_{0.92f};

    // adaptive notch filter
    dsp::AdaptiveNotchFilter<float> adaptiveNotch_{0.996f, 0.998f};
    decltype(adaptiveNotch_)::State adaptiveNotchState_{};
    dsp::NotchFilter<float> notch_{0.98f};
    decltype(notch_)::State notchState_{};
    dsp::OnePole<float> ampFollower_{};
    decltype(ampFollower_)::State ampFollowerState_{};
    decltype(ampFollower_)::State poleFollowerState_{};

    bool voiced_{};
    int nConfirmPitch_{};

    dsp::Buffer<float, nextTo(nsdfDL_) + 512> bufferDec_;
    std::array<float, decltype(bufferDec_)::kSize> bufferDecData_{};
};
} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
