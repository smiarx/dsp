#pragma once

#include "../../Buffer.h"
#include "../../Delay.h"
#include "../../LFO.h"
#include "../../TapeDelay.h"

static constexpr auto N = 2;

class TapeDelay
{
  public:
    static constexpr auto MaxBlockSize = 512;
    static constexpr int MaxDelay      = 48000 * 1.0;

    static constexpr auto speedSmoothTime = 1.0f;
    static constexpr auto speedModFreq    = 0.242f;

    TapeDelay(float sampleRate)
    {
        buffer_.setBuffer(bufferArray_);
        setSampleRate(sampleRate);

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

    void update(float delay, float feedback, float drywet);
    void process(float **__restrict in, float **__restrict out, int count);

  private:
    float sampleRate_{48000.f};
    float freqScale_{2.f / 48000.f};
    float invSampleRate_{1.f / 48000.f};

    float delay_{0.f};
    float feedback_{0.f};
    float drywet_{1.f};

    using TapePosition = dsp::TapePosition<MaxDelay>;
    float targetSpeed_{0};
    float speed_{0};
    float speedSmooth_{0.f};
    float speedMod_{0.000004f};
    dsp::LFOParabolic<1> speedLFO_;
    TapePosition tapePos_;
    dsp::TapTape tapTape_;
    dsp::DelayLine<MaxDelay> delayline_;

    static constexpr auto BufferSize = nextTo(delayline_);
    dsp::Buffer<dsp::Signal<2>, BufferSize> buffer_;
    dsp::Signal<N> bufferArray_[decltype(buffer_)::Size] = {{0.f}};

    // write tmp buffers
    dsp::Signal<N> x_[MaxBlockSize] = {{0.f}};
};
