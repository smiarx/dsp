#pragma once

#include <dsp/control/Smoother.hpp>
#include <dsp/core/Buffer.hpp>
#include <dsp/core/SampleRate.hpp>
#include <dsp/delay/Delay.hpp>
#include <dsp/delay/TapeDelay.hpp>
#include <dsp/filters/VAFilters.hpp>
#include <dsp/modulation/LFO.hpp>
#include <dsp/resampling/Kernel.hpp>
#include <dsp/resampling/Windows.hpp>

#ifdef TAPEDELAY_SWITCH_INDICATOR
#include <atomic>
#endif

#include "TapeDelayDefines.h"

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

using namespace dsp::units;

class TapeDelay
{
  public:
    static constexpr auto kN = 2;
    using type               = float;
    using mtype              = dsp::multi<type, kN>;

    using Hz      = dsp::freq::Hz<type>;
    using seconds = dsp::time::Seconds<type>;

    static constexpr auto kMaxBlockSize      = 512;
    static constexpr auto kDefaultSampleRate = dsp::SampleRate{48000._Hzf};
    static constexpr auto kMaxDelay =
        dsp::time::Seconds<type>(TAPEDELAY_MAX_DELAY);

    static constexpr type kReverseDelayMaxRatio =
        type(TAPEDELAY_REVERSE_MAX_RATIO);

    static constexpr auto kSpeedSmoothTime = 0.7_sf;
    static constexpr auto kSpeedModFreq    = 0.242_Hzf;
    static constexpr auto kSpeedModAmp     = 0.02f;

    static constexpr auto kKernelSize = 4;

    static constexpr auto kDelayBufSize =
        static_cast<size_t>(kDefaultSampleRate.secsToSamples(kMaxDelay).get());

    // lookup table for cross fading between two taps
    static constexpr auto kFadeSize = 2048;
    class FadeLut : public std::array<float, kFadeSize>
    {
      public:
        FadeLut()
        {
            for (size_t n = 0; n < kFadeSize; ++n) {
                float x =
                    (static_cast<float>(n) + 1.f) / (float(kFadeSize) + 1.f);
                (*this)[n] = dsp::windows::Hann::generate(x);
            }
        }
    };
    // we use pointer instead of variable so that it's not initialized at
    // startup - crashes with bad instruction arch if avx is not available
    static FadeLut *fadeLut;

    enum Mode {
        kNormal    = 0,
        kBackForth = 1,
        kReverse   = 2,
    };

    TapeDelay() = default;

    template <class ReAlloc = decltype(std::realloc)>
    void prepare(float sampleRate, int blockSize,
                 ReAlloc realloc = std::realloc);

    template <class Free = decltype(std::free)>
    void free(Free free = std::free);

    // getters
    [[nodiscard]] float getDelay() const { return delay_.get(); }
    [[nodiscard]] float getDrift() const { return drift_; }
    [[nodiscard]] float getCutLowPass() const { return cutlowpass_.get(); }
    [[nodiscard]] float getCutHiPass() const { return cuthighpass_.get(); }
    [[nodiscard]] float getSaturation() const
    {
        return dsp::getlane<0>(saturation_.getTarget());
    }
    [[nodiscard]] float getFeedback() const { return feedback_; }
    [[nodiscard]] float getDryWet() const { return drywet_; }
    [[nodiscard]] Mode getMode() const { return mode_; }

    // setters
    void setDelay(seconds delay, int blockSize);
    void setDrift(float drift, int blockSize);
    void setCutLowPass(Hz cutlowpass, int blockSize);
    void setCutHiPass(Hz cuthighpass, int blockSize);
    void setSaturation(float saturation, int blockSize);
    void setFeedback(float feedback, int blockSize);
    void setDryWet(float drywet, int blockSize);
    void setMode(Mode mode, int blockSize);

    void update(seconds delay, float feedback, Hz cutlowpass, Hz cuthighpass,
                float saturation, float drift, Mode mode, float drywet,
                int blockSize);
    void process(const float *const *__restrict in,
                 float *const *__restrict out, int count);

  private:
    dsp::SampleRate<type> sampleRate_{kDefaultSampleRate};
    int maxBlockSize_{};
    int maxBlockSizeWithDelay_{};

    seconds delay_{0._sf};
    type feedback_{};
    dsp::ControlSmoother<mtype, true> feedbackCompensated_{};
    type drywet_{};
    dsp::ControlSmoother<mtype, true> dry_{1.f};
    dsp::ControlSmoother<mtype, true> wet_{};

    // tape movement
    using TapePosition = dsp::TapePosition<kDelayBufSize>;
    type targetSpeed_{};
    type speed_{};
    type speedSmooth_{};

    TapePosition tapePos_;
    using TapeInterp = dsp::TapKernel<
        mtype, dsp::kernels::Sinc<kKernelSize, dsp::windows::Kaiser<140>>, 64>;
    // using TapeInterp = dsp::TapNoInterp<mtype>;
    using TapTape = dsp::TapTape<TapeInterp>;
    TapTape tapTape_[2];
    TapTape tapReverse_;
    TapePosition::position_t reverseDist_[2]{0, 0};
    int tapId_{0};
    int fadePos_{-1};

    // switch tape width mode
    void switchTap(Mode mode);
#ifdef TAPEDELAY_SWITCH_INDICATOR
  public:
    std::atomic<bool> &getSwitchIndicator() { return switchIndicator_; }

  private:
    std::atomic<bool> switchIndicator_;
#endif
    // read function
    template <Mode, bool check, class Ctxt>
    bool read(Ctxt ctxt, int tapId, TapePosition::position_t speed);
    template <Mode, class Ctxt> int readBlock(Ctxt ctxt);
    // move tape function
    template <class Ctxt> TapePosition::position_t moveTape(Ctxt);

    dsp::DelayLine<kDelayBufSize / 3> delaylineReverse_;
    dsp::DelayLine<kDelayBufSize * 2 / 3, nextTo(delaylineReverse_)> delayline_;

    // mode
    Mode mode_{kNormal};
    Mode oldMode_{kNormal};

    // speed modulation
    float drift_{0.f};
    dsp::ControlSmoother<type> speedMod_{};
    dsp::lfo::Sine<type> speedLFO_;

    // saturation
    dsp::ControlSmoother<mtype, true> saturation_{};
    mtype pregain_{1};
    mtype postgain_{1};

    // low pass filter
    Hz cutlowpass_{0._Hzf};
    dsp::va::SVF<mtype, dsp::va::kLowPass> lpf_{};
    decltype(lpf_)::State lpfMem_{};

    // high pass filter
    Hz cuthighpass_{0._Hzf};
    dsp::va::SVF<mtype, dsp::va::kHighPass> hpf_{};
    decltype(hpf_)::State hpfMem_{};

    static constexpr auto kBufferSize = nextTo(delayline_);
    dsp::Buffer<mtype, kBufferSize> buffer_;

    // tmp buffers
    mtype *__restrict x_{nullptr};
};

template <class ReAlloc>
void TapeDelay::prepare(float sampleRate, int blockSize, ReAlloc realloc)
{
    sampleRate_   = dsp::SampleRate{dsp::freq::Hz{sampleRate}};
    maxBlockSize_ = std::min(blockSize, kMaxBlockSize);

    speedSmooth_ =
        1.f - std::pow(0.001f,
                       1.f / sampleRate_.secsToSamples(kSpeedSmoothTime).get());

    speedLFO_.setFreq(sampleRate_.hzToNyquistRatio(kSpeedModFreq));

    // alloc ressources
    x_ = (mtype *)realloc(x_,
                          sizeof(mtype) * static_cast<size_t>(maxBlockSize_));

    auto *buf                      = buffer_.getData();
    constexpr auto kRealBufferSize = sizeof(mtype) * decltype(buffer_)::kSize;
    buf                            = (mtype *)realloc(buf, kRealBufferSize);
    memset(buf, 0, kRealBufferSize);
    buffer_.setData(buf);

    if (fadeLut == nullptr) {
        fadeLut = (FadeLut *)realloc(fadeLut, sizeof(FadeLut));
        new (fadeLut) FadeLut{};
    }
}

template <class Free> void TapeDelay::free(Free free)
{
    free(x_);
    x_ = nullptr;

    auto *buf = buffer_.getData();
    free(buf);
    buffer_.setData(nullptr);
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
