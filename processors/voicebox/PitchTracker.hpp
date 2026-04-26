#pragma once

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

  public:
    PitchTracker()
    {
        buffer_.setData(bufferData_.data());
        bufferDec_.setData(bufferDecData_.data());
    };

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

    // amp follower
    dsp::OnePole<float> ampFollower_{};
    decltype(ampFollower_)::State ampFollowerState_{};

    // nsdf
    dsp::NSDF<float, 120> nsdf_{};
    decltype(nsdf_)::DL<0> nsdfDL_{};
    float downPeriod_{};
    float period_{};

    enum HarmonicJump {
        kJumpOctaveDown,
        kNoJump,
        kJumpOctaveUp,
        kJump3rdHarmonicUp,
        kNumJumps,
    } harmonicJumpType_{kNoJump};
    int harmonicJumpCounter_{};

    // zero crossing voicing detection
    dsp::ZeroCrossing<float> zeroCrossing_{};
    dsp::OnePole<float> zeroCrossingFollower_{};
    decltype(zeroCrossingFollower_)::State zeroCrossingFollowerState_{};
    float zeroCrossingBuffer_[kMaxBlockSize]{};
    // voicing detection
    dsp::HMM2<float> hmm2_{0.92f};

    // track
    dsp::NSDFTrack<float, 32, 960> nsdfTrack_{};
    decltype(nsdfTrack_)::DL<0> nsdfTrackDL_{};

    bool voiced_{};
    int nConfirmPitch_{};

    dsp::Buffer<float, nextTo(nsdfDL_) + kMaxBlockSize> bufferDec_;
    std::array<float, decltype(bufferDec_)::kSize> bufferDecData_{};

    dsp::Buffer<float, nextTo(nsdfTrackDL_) + kMaxBlockSize> buffer_;
    std::array<float, decltype(buffer_)::kSize> bufferData_{};
};

// class PitchTracker
//{
//     static constexpr auto kMaxBlockSize       = 512;
//     static constexpr auto kDecFactor = 8;
//
//   public:
//     PitchTracker() { bufferDec_.setData(bufferDecData_.data()); };
//
//     void process(const float *const *__restrict in,
//                  float *const *__restrict out, int count);
//
//     void prepare(float sampleRate);
//
//     float findNSDFPeriod(float zeroCrossing);
//
//   private:
//     dsp::va::OnePole<float, dsp::va::kLowPass> softLowpass_{};
//     decltype(softLowpass_)::State softLowpassState_{};
//
//     // decimate filter
//     dsp::IIRFilter<float, 8> nyquistFilter_{};
//     decltype(nyquistFilter_)::State nyquistFilterState_{};
//
//     // decimate
//     int decimateId_{};
//
//     // nsdf
//     dsp::NSDF<float, 120> nsdf_{};
//     decltype(nsdf_)::DL<0> nsdfDL_{};
//
//     // zero crossing voicing detection
//     dsp::ZeroCrossing<float> zeroCrossing_{};
//     dsp::OnePole<float> zeroCrossingFollower_{};
//     decltype(zeroCrossingFollower_)::State zeroCrossingFollowerState_{};
//     float zeroCrossingBuffer_[kMaxBlockSize]{};
//     // voicing detection
//     dsp::HMM2<float> hmm2_{0.92f};
//
//     // adaptive notch filter
//     dsp::AdaptiveNotchFilter<float> adaptiveNotch_{0.990f, 0.9920f};
//     decltype(adaptiveNotch_)::State adaptiveNotchState_{};
//     float lp_[2];
//     float lp2_[2];
//     dsp::NotchFilter<float> notch_{0.99f};
//     decltype(notch_)::State notchState_{};
//     decltype(notch_)::State notchState2_{};
//     dsp::OnePole<float> ampFollower_{};
//     decltype(ampFollower_)::State ampFollowerState_{};
//     decltype(ampFollower_)::State poleFollowerState_{};
//
//     bool voiced_{};
//     int nConfirmPitch_{};
//
//     dsp::Buffer<float, nextTo(nsdfDL_) + 512> bufferDec_;
//     std::array<float, decltype(bufferDec_)::kSize> bufferDecData_{};
// };

// class PitchTracker
//{
//     static constexpr auto kMaxBlockSize       = 512;
//
//   public:
//     PitchTracker() = default;
//
//     void process(const float *const *__restrict in,
//                  float *const *__restrict out, int count);
//
//     void prepare(float sampleRate);
//
//   private:
//     dsp::va::SVF<float, dsp::va::kLowPass> softLowpass_{};
//     decltype(softLowpass_)::State softLowpassState_{};
//
//     // decimate filter
//     dsp::IIRFilter<float, 8> nyquistFilter_{};
//     decltype(nyquistFilter_)::State nyquistFilterState_{};
//
//     // adaptive notch filter
//     dsp::AdaptiveNotchFilter<float> adaptiveNotch_{0.98f, 0.990f};
//     decltype(adaptiveNotch_)::State adaptiveNotchState_{};
//     decltype(adaptiveNotch_)::State adaptiveNotchState2_{};
//     dsp::NotchFilter<float> notch_{0.98f};
//     decltype(notch_)::State notchState_{};
//
//     float lp_[2];
//
//     dsp::IIRFilter<float,4> dcRemove_{};
//     decltype(dcRemove_)::State dcRemoveState_{};
//
//     // zero crossing voicing detection
//     dsp::ZeroCrossing<float> zeroCrossing_{};
//     dsp::OnePole<float> zeroCrossingFollower_{};
//     decltype(zeroCrossingFollower_)::State zeroCrossingFollowerState_{};
//     float zeroCrossingBuffer_[kMaxBlockSize]{};
//     // voicing detection
//     dsp::HMM2<float> hmm2_{0.99f};
//     float prevA_{1.99f};
// };
} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
