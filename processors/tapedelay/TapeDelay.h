#pragma once

#include "dsp/Buffer.h"
#include "dsp/Delay.h"
#include "dsp/Kernel.h"
#include "dsp/LFO.h"
#include "dsp/Smoother.h"
#include "dsp/TapeDelay.h"
#include "dsp/VAFilters.h"
#include "dsp/Window.h"

#ifdef TAPEDELAY_SWITCH_INDICATOR
#include <atomic>
#endif

namespace processors
{

class TapeDelay
{
  public:
    static constexpr auto kN = 2;

    static constexpr auto kMaxBlockSize      = 512;
    static constexpr auto kDefaultSampleRate = 48000.f;
    static constexpr int kMaxDelay           = 5.f;

    static constexpr float kReverseDelayMaxRatio = 3.1;

    static constexpr auto kSpeedSmoothTime = 0.7f;
    static constexpr auto kSpeedModFreq    = 0.242f;
    static constexpr auto kSpeedModAmp     = 0.013f;

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
                (*this)[n] = dsp::window::Hann::generate(x);
            }
        }
    };
    static FadeLut fadeLut;

    enum Mode {
        kNormal    = 0,
        kBackForth = 1,
        kReverse   = 2,
    };

    TapeDelay()
    {
        inFor(pregain_, k, i)
        {
            pregain_[k][i]  = 1.f;
            postgain_[k][i] = 1.f;
        }
    }

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
        return saturation_.getTarget()[0];
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
    float freqScale_{2.f / kDefaultSampleRate};
    float sampleRate_{1.f / kDefaultSampleRate};
    float invSampleRate_{1.f / kDefaultSampleRate};
    int maxBlockSize_{};

    float delay_{0.f};
    float feedback_{0.f};
    dsp::ControlSmoother<2, true> feedbackCompensated_{{0.f, 0.f}};
    float drywet_{0.f};
    dsp::ControlSmoother<2> dry_{{1.f, 1.f}};
    dsp::ControlSmoother<2> wet_{{0.f, 0.f}};

    // tape movement
    using TapePosition = dsp::TapePosition<kDelayBufSize>;
    float targetSpeed_{0};
    float speed_{0};
    float speedSmooth_{0.f};

    TapePosition tapePos_;
    using TapeInterp = dsp::TapKernel<
        2, dsp::kernel::Sinc<kKernelSize, dsp::window::Kaiser<140>>, 64>;
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
    TapePosition::position_t moveTape();

    dsp::DelayLine<kDelayBufSize / 3> delaylineReverse_;
    dsp::DelayLine<kDelayBufSize * 2 / 3, nextTo(delaylineReverse_)> delayline_;

    // mode
    Mode mode_{kNormal};
    Mode oldMode_{kNormal};

    // speed modulation
    float drift_{0.f};
    dsp::ControlSmoother<1> speedMod_{{0.f}};
    dsp::LFOSine<1> speedLFO_;

    // saturation
    dsp::ControlSmoother<1> saturation_{{0.f}};
    dsp::fSample<2>::Vector pregain_;
    dsp::fSample<2>::Vector postgain_;

    // low pass filter
    float cutlowpass_{0.f};
    dsp::va::SVF<kN, dsp::va::kLowPass> lpf_{};
    decltype(lpf_)::State lpfMem_{};

    // high pass filter
    float cuthighpass_{0.f};
    dsp::va::SVF<kN, dsp::va::kHighPass> hpf_{};
    decltype(hpf_)::State hpfMem_{};

    static constexpr auto kBufferSize = nextTo(delayline_);
    dsp::Buffer<dsp::fSample<kN>, kBufferSize> buffer_;

    // tmp buffers
    dsp::fSample<kN> *__restrict x_{nullptr};
};

template <class ReAlloc>
void TapeDelay::prepare(float sampleRate, int blockSize, ReAlloc realloc)
{
    sampleRate_    = sampleRate;
    invSampleRate_ = 1.f / sampleRate;
    freqScale_     = 2.f * invSampleRate_;
    maxBlockSize_  = std::min(blockSize, kMaxBlockSize);

    speedSmooth_ =
        1.f - std::pow(0.001f, 1.f / kSpeedSmoothTime * invSampleRate_);

    speedLFO_.setFreq({freqScale_ * kSpeedModFreq});

    // alloc ressources
    x_ = (dsp::fSample<kN> *)realloc(
        x_, sizeof(dsp::fSample<kN>) * static_cast<size_t>(maxBlockSize_));

    auto *buf = buffer_.getBuffer();
    constexpr auto kRealBufferSize =
        sizeof(dsp::fSample<kN>) * decltype(buffer_)::kSize;
    buf = (dsp::fSample<kN> *)realloc(buf, kRealBufferSize);
    memset(buf, 0, kRealBufferSize);
    buffer_.setBuffer(buf);
}

template <class Free> void TapeDelay::free(Free free)
{
    free(x_);
    x_ = nullptr;

    auto *buf = buffer_.getBuffer();
    free(buf);
    buffer_.setBuffer(nullptr);
}

} // namespace processors
