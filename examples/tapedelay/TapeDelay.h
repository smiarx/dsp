#pragma once

#include "../../Buffer.h"
#include "../../Delay.h"
#include "../../TapeDelay.h"

static constexpr auto N = 2;

class TapeDelay
{
  public:
    static constexpr auto MaxBlockSize = 512;
    static constexpr int MaxDelay      = 48000 * 1.0;

    TapeDelay(float sampleRate)
    {
        buffer_.setBuffer(bufferArray_);
        setSampleRate(sampleRate);
    }

    void setSampleRate(float sR)
    {
        sampleRate_    = sR;
        invSampleRate_ = 1.f / sR;
        freqScale_     = 2.f * invSampleRate_;
    }

    void update(float delay);
    void process(float **__restrict in, float **__restrict out, int count);

  private:
    float sampleRate_{48000.f};
    float freqScale_{2.f / 48000.f};
    float invSampleRate_{1.f / 48000.f};

    float delay_{0.f};

    using TapePosition = dsp::TapePosition<MaxDelay>;
    TapePosition::position_t targetSpeed_{0};
    TapePosition::position_t speed_{0};
    TapePosition tapePos_;
    dsp::TapTape tapTape_;
    dsp::DelayLine<MaxDelay> delayline_;

    static constexpr auto BufferSize = nextTo(delayline_);
    dsp::Buffer<dsp::Signal<2>, BufferSize> buffer_;
    dsp::Signal<N> bufferArray_[decltype(buffer_)::Size] = {{0.f}};
};
