#pragma once

#include "dsp/Buffer.h"
#include "dsp/Delay.h"
#include "dsp/Kernel.h"
#include "dsp/LFO.h"
#include "dsp/Smoother.h"
#include "dsp/TapeDelay.h"
#include "dsp/VAFilters.h"
#include "dsp/Windows.h"

#ifdef TAPEDELAY_SWITCH_INDICATOR
#include <atomic>
#endif

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

class TapeDelay
{
  public:
    static constexpr auto kN = 2;
    using type               = float;
    using mtype              = dsp::multi<type, kN>;

    static constexpr auto kMaxBlockSize      = 512;
    static constexpr type kDefaultSampleRate = 48000;
    static constexpr type kMaxDelay          = 5;

    static constexpr type kReverseDelayMaxRatio = 3.1;

    static constexpr type kSpeedSmoothTime = 0.7;
    static constexpr type kSpeedModFreq    = 0.242;
    static constexpr type kSpeedModAmp     = 0.02;

    static constexpr auto kKernelSize = 4;

    static constexpr auto kDelayBufSize =
        static_cast<size_t>(kDefaultSampleRate * kMaxDelay);

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
    static FadeLut fadeLut;

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
    [[nodiscard]] float getDelay() const { return delay_; }
    [[nodiscard]] float getDrift() const { return drift_; }
    [[nodiscard]] float getCutLowPass() const { return cutlowpass_; }
    [[nodiscard]] float getCutHiPass() const { return cuthighpass_; }
    [[nodiscard]] float getSaturation() const
    {
        return dsp::getlane<0>(saturation_.getTarget());
    }
    [[nodiscard]] float getFeedback() const { return feedback_; }
    [[nodiscard]] float getDryWet() const { return drywet_; }
    [[nodiscard]] Mode getMode() const { return mode_; }

    // setters
    void setDelay(float delay, int blockSize);
    void setDrift(float drift, int blockSize);
    void setCutLowPass(float cutlowpass, int blockSize);
    void setCutHiPass(float cuthighpass, int blockSize);
    void setSaturation(float saturation, int blockSize);
    void setFeedback(float feedback, int blockSize);
    void setDryWet(float drywet, int blockSize);
    void setMode(Mode mode, int blockSize);

    void update(float delay, float feedback, float cutlowpass,
                float cuthighpass, float saturation, float drift, Mode mode,
                float drywet, int blockSize);
    void process(const float *const *__restrict in,
                 float *const *__restrict out, int count);

  private:
    type freqScale_{2. / kDefaultSampleRate};
    type sampleRate_{1. / kDefaultSampleRate};
    int maxBlockSize_{};
    int maxBlockSizeWithDelay_{};

    type delay_{};
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
    type cutlowpass_{};
    dsp::va::SVF<mtype, dsp::va::kLowPass> lpf_{};
    decltype(lpf_)::State lpfMem_{};

    // high pass filter
    float cuthighpass_{0.f};
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
    sampleRate_        = sampleRate;
    auto invSampleRate = 1.f / sampleRate;
    freqScale_         = 2.f * invSampleRate;
    maxBlockSize_      = std::min(blockSize, kMaxBlockSize);

    speedSmooth_ =
        1.f - std::pow(0.001f, 1.f / kSpeedSmoothTime * invSampleRate);

    speedLFO_.setFreq({freqScale_ * kSpeedModFreq});

    // alloc ressources
    x_ = (mtype *)realloc(x_,
                          sizeof(mtype) * static_cast<size_t>(maxBlockSize_));

    auto *buf                      = buffer_.getData();
    constexpr auto kRealBufferSize = sizeof(mtype) * decltype(buffer_)::kSize;
    buf                            = (mtype *)realloc(buf, kRealBufferSize);
    memset(buf, 0, kRealBufferSize);
    buffer_.setData(buf);
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
