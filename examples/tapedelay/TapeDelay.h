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
        sampleRate_ = sR;
        freqScale_  = 2.f / sR;
    }

    void update(float delay);
    void process(float **__restrict in, float **__restrict out, int count);

  private:
    float sampleRate_{48000.f};
    float freqScale_{2.f / 48000.f};

    float speed_{1.f / 0.5f / 48000.f};

    dsp::DelayLine<MaxDelay> delayline_;
    dsp::TapePosition<MaxDelay> tapePos_;
    dsp::TapTape tapTape_;

    static constexpr auto BufferSize = nextTo(delayline_);
    dsp::Buffer<dsp::Signal<2>, BufferSize> buffer_;
    dsp::Signal<N> bufferArray_[decltype(buffer_)::Size] = {{0.f}};
};
