#pragma once

#include "dsp/Buffer.h"
#include "dsp/Delay.h"
#include "dsp/Kernel.h"
#include "dsp/LFO.h"
#include "dsp/Smoother.h"
#include "dsp/TapeDelay.h"
#include "dsp/VAFilters.h"
#include "dsp/Window.h"

namespace processors
{

class TapeDelay
{
  public:
    static constexpr auto N = 2;

    static constexpr auto MaxBlockSize      = 512;
    static constexpr auto DefaultSampleRate = 48000.f;
    static constexpr int MaxDelay           = 2.f;

    static constexpr auto speedSmoothTime = 0.7f;
    static constexpr auto speedModFreq    = 0.242f;
    static constexpr auto speedModAmp     = 0.013f;

    static constexpr auto KernelSize = 4;

    static constexpr size_t DelayBufSize = DefaultSampleRate * MaxDelay;

    // lookup table for cross fading between two taps
    static constexpr auto FadeSize = 2048;
    class FadeLut : public std::array<float, FadeSize>
    {
      public:
        FadeLut()
        {
            for (size_t n = 0; n < FadeSize; ++n) {
                (*this)[n] =
                    dsp::window::Hann::generate((n + 1.f) / (FadeSize + 1.f));
            }
        }
    };
    static FadeLut fadeLut;

    enum Mode {
        Normal    = 0,
        BackForth = 1,
        Reverse   = 2,
    };

    TapeDelay()
    {
        buffer_.setBuffer(bufferArray_);

        inFor(pregain_, k, i)
        {
            pregain_[k][i]  = 1.f;
            postgain_[k][i] = 1.f;
        }
    }

    void setSampleRate(float sR)
    {
        sampleRate_    = sR;
        invSampleRate_ = 1.f / sR;
        freqScale_     = 2.f * invSampleRate_;

        speedSmooth_ =
            1.f - powf(0.001, 1.f / speedSmoothTime * invSampleRate_);

        speedLFO_.setFreq({freqScale_ * speedModFreq});
    }

    void setBlockSize(int blockSize)
    {
        blockSize_    = blockSize;
        invBlockSize_ = 1.f / blockSize;
        blockRate_    = sampleRate_ * invBlockSize_;
    }

    // getters
    float getDelay() const { return delay_; }
    float getDrift() const { return drift_; }
    float getCutLowPass() const { return cutlowpass_; }
    float getCutHiPass() const { return cuthighpass_; }
    float getSaturation() const { return saturation_.getTarget()[0]; }
    float getFeedback() const { return feedback_; }
    float getDryWet() const { return drywet_.getTarget()[0]; }

    // setters
    void setDelay(float delay);
    void setDrift(float drift);
    void setCutLowPass(float cutlowpass);
    void setCutHiPass(float cuthipass);
    void setSaturation(float saturation);
    void setFeedback(float feedback);
    void setDryWet(float drywet);
    void setMode(Mode mode);

    void update(float delay, float feedback, float cutlp, float cuthp,
                float saturation, float drift, Mode mode, float drywet);
    void process(const float *const *__restrict in,
                 float *const *__restrict out, int count);

  private:
    float sampleRate_{48000.f};
    float freqScale_{2.f / 48000.f};
    float invSampleRate_{1.f / 48000.f};
    int blockSize_{0};
    float invBlockSize_{0.f};
    float blockRate_{0.f};

    float delay_{0.f};
    float feedback_{0.f};
    dsp::ControlSmoother<2, true> feedbackCompensated_{{0.f, 0.f}};
    dsp::ControlSmoother<2, true> drywet_{{0.f, 0.f}};

    // tape movement
    using TapePosition = dsp::TapePosition<DelayBufSize>;
    float targetSpeed_{0};
    float speed_{0};
    float speedSmooth_{0.f};

    TapePosition tapePos_;
    using TapeInterp = dsp::TapKernel<
        2, dsp::kernel::Sinc<KernelSize, dsp::window::Kaiser<140>>, 64>;
    using TapTape = dsp::TapTape<TapeInterp>;
    TapTape tapTape_[2];
    TapTape tapReverse_;
    TapePosition::position_t reverseDist_[2]{0, 0};
    size_t tapId_{0};
    int fadePos_{-1};

    // switch tape width mode
    void switchTap(Mode mode);
    // read function
    template <Mode, bool check, class Ctxt>
    bool read(Ctxt ctxt, int tapId, TapePosition::position_t speed);
    template <Mode, class Ctxt> int readBlock(Ctxt ctxt);
    // move tape function
    float moveTape();

    dsp::DelayLine<DelayBufSize / 3> delaylineReverse_;
    dsp::DelayLine<DelayBufSize * 2 / 3, nextTo(delaylineReverse_)> delayline_;

    // mode
    Mode mode_{Normal};
    Mode oldMode_{Normal};

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
    dsp::va::SVF<N, dsp::va::LowPass> lpf_{};
    decltype(lpf_)::State lpfMem_{};

    // high pass filter
    float cuthighpass_{0.f};
    dsp::va::SVF<N, dsp::va::HighPass> hpf_{};
    decltype(hpf_)::State hpfMem_{};

    static constexpr auto BufferSize = nextTo(delayline_);
    dsp::Buffer<dsp::fSample<2>, BufferSize> buffer_;
    dsp::fSample<N> bufferArray_[decltype(buffer_)::Size] = {{0.f}};

    // write tmp buffers
    dsp::fSample<N> x_[MaxBlockSize] = {{0.f}};
};
} // namespace processors
