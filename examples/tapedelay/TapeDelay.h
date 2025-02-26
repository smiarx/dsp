#pragma once

#include "../../Buffer.h"
#include "../../Delay.h"
#include "../../Kernel.h"
#include "../../LFO.h"
#include "../../Smoother.h"
#include "../../TapeDelay.h"
#include "../../VAFilters.h"
#include "../../Window.h"

static constexpr auto N = 2;

class TapeDelay
{
  public:
    static constexpr auto MaxBlockSize = 512;
    static constexpr int MaxDelay      = 48000 * 1.0;

    static constexpr auto speedSmoothTime = 0.7f;
    static constexpr auto speedModFreq    = 0.242f;
    static constexpr auto speedModAmp     = 0.013f;

    static constexpr auto KernelSize = 4;

    // lookup table for cross fading between two taps
    static constexpr auto FadeSize = 8;
    class FadeLut : public std::array<float, FadeSize>
    {
      public:
        FadeLut()
        {
            for (size_t n = 0; n < size(); ++n) {
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

    TapeDelay(float sampleRate, int blockSize)
    {
        buffer_.setBuffer(bufferArray_);
        setSampleRate(sampleRate);
        setBlockSize(blockSize);

        inFor(pregain_, k, i)
        {
            pregain_[k][i]  = 1.f;
            postgain_[k][i] = 1.f;
        }

        speedLFO_.setFreq({freqScale_ * speedModFreq});
    }

    void setSampleRate(float sR)
    {
        sampleRate_    = sR;
        invSampleRate_ = 1.f / sR;
        freqScale_     = 2.f * invSampleRate_;

        speedSmooth_ =
            1.f - powf(0.001, 1.f / speedSmoothTime * invSampleRate_);
    }

    void setBlockSize(int blockSize)
    {
        blockSize_    = blockSize;
        invBlockSize_ = 1.f / blockSize;
        blockRate_    = sampleRate_ * invBlockSize_;
    }

    void update(float delay, float feedback, float cutlp, float cuthp,
                float saturation, float drift, Mode mode, float drywet);
    void process(float **__restrict in, float **__restrict out, int count);

  private:
    float sampleRate_{48000.f};
    float freqScale_{2.f / 48000.f};
    float invSampleRate_{1.f / 48000.f};
    int blockSize_{0};
    float invBlockSize_{0.f};
    float blockRate_{0.f};

    float delay_{0.f};
    dsp::ControlSmoother<2, true> feedback_{{0.f, 0.f}};
    dsp::ControlSmoother<2, true> drywet_{{0.f, 0.f}};

    // tape movement
    using TapePosition = dsp::TapePosition<MaxDelay>;
    float targetSpeed_{0};
    float speed_{0};
    float speedSmooth_{0.f};

    TapePosition tapePos_;
    using TapeInterp = dsp::TapKernel<
        2, dsp::kernel::Sinc<KernelSize, dsp::window::Kaiser<140>>, 64>;
    using TapTape = dsp::TapTape<TapeInterp>;
    TapTape tapTape_[2];
    TapePosition::position_t reverseDist_[2]{0, 0};
    size_t tapId_{0};
    int fadePos_{-1};

    // read function
    template <Mode, bool check, class Ctxt>
    bool read(Ctxt ctxt, int tapId, TapePosition::position_t speed);
    template <Mode, class Ctxt> int readBlock(Ctxt ctxt);
    void switchTap(Mode mode);

    dsp::DelayLine<MaxDelay> delayline_;

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
    dsp::va::SVF<N, dsp::va::LowPass> lpf_;
    decltype(lpf_)::State lpfMem_;

    // high pass filter
    float cuthighpass_{0.f};
    dsp::va::SVF<N, dsp::va::HighPass> hpf_;
    decltype(hpf_)::State hpfMem_;

    static constexpr auto BufferSize = nextTo(delayline_);
    dsp::Buffer<dsp::fSample<2>, BufferSize> buffer_;
    dsp::fSample<N> bufferArray_[decltype(buffer_)::Size] = {{0.f}};

    // write tmp buffers
    dsp::fSample<N> x_[MaxBlockSize] = {{0.f}};
};
