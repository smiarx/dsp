

#include "../../AllPass.h"
#include "../../Filter.h"
#include "../../LFO.h"
#include "../../MultiRate.h"
#include "../../Noise.h"

static constexpr auto N = 4;

class Springs
{
  public:
    static constexpr auto MaxBlockSize = 512;
    static constexpr auto MaxDecimate  = 8;
    static constexpr auto CascadeL     = 80;

    static constexpr float freqFactor[] = {0.97743f,1.04391f,1.0593f,0.934f};
    static constexpr float RFactor[] = {1.07743f,0.94391f,0.9893f,1.034f};
    static constexpr float loopTdFactor[] = {0.97874f,1.03913f,0.953872f,1.18373f};
    static constexpr float loopModFreq[] = {0.35f,0.564f,0.46,0.20f};
    static constexpr float loopModFactor[] = {0.002743f,0.00205f,0.00348f,0.0021f};
    static constexpr float loopEchoGain = 0.076f;
    static constexpr float loopRippleGain = 0.036f;

    Springs(float sampleRate)
    {
        buffer_.setBuffer(bufferArray_);
        bufferDec_.setBuffer(bufferDecArray_);
        setSampleRate(sampleRate);

        /* loop modulation */
        dsp::iNoise<N> noise;
        loopMod_.setPhase(noise.process());

        dsp::Signal<N>freq;
        for(int i = 0; i < N; ++i){
            freq[i] = loopModFreq[i] * freqScale_;
        }
        loopMod_.setFreq(freq);
    }

    int M_{1};
    float sampleRate_{48000.f};
    float freqScale_{2.f / 48000.f};
    float R_{0.f};
    float freq_{1.f};
    float Td_{0.f};
    float T60_{0.f};

    dsp::AllPass2<N> allpass_;
    typename decltype(allpass_)::DL allpassdl_[CascadeL];

    dsp::IIRFilter<N, 4> lowpass_;
    typename decltype(lowpass_)::DL<N> lowpassdl_;

    using MR = dsp::MultiRate<N, 15, MaxDecimate>;
    MR::DLDecimate dldecimate_;
    MR::DLInterpolate dlinterpolate_;
    int decimateId_{0};

    static constexpr auto BufSize = nextTo(dldecimate_) + MaxBlockSize;

    using MRs = MR::WithBuffer<BufSize>;
    const MRs::Base *multirate_;

    dsp::Buffer<dsp::Signal<N>, BufSize> buffer_;
    dsp::Signal<N> bufferArray_[decltype(buffer_)::Size] = {{0.f}};

    static constexpr int LoopLength = 0.25 * 48000;
    float loopGain_{0.f};
    dsp::Signal<N> loopTd_{0.f};
    dsp::Signal<N> loopModAmp_{0.f};
    dsp::LFOParabolic<N> loopMod_;
    dsp::TapNoInterp<N> loopEcho_;
    using LoopType = dsp::TapAllPass<N, dsp::TapNoInterp<N>>;
    dsp::NestedDelayLine<dsp::DelayLine<LoopLength>, dsp::CopyDelayLine<N, 1>>
        loopdl_;
    dsp::DelayLine<LoopLength/5, nextTo(loopdl_)> loopEchoDL_;
    dsp::CopyDelayLine<N,1,nextTo(loopEchoDL_)> loopRippleDL_;

    static constexpr auto BufDecSize = nextTo(loopRippleDL_) + MaxBlockSize;
    dsp::Buffer<dsp::Signal<N>, BufDecSize> bufferDec_;
    dsp::Signal<N> bufferDecArray_[decltype(bufferDec_)::Size] = {{0.f}};

    dsp::Signal<N> x_[MaxBlockSize];
    dsp::Signal<N> xdecimate_[MaxBlockSize / 2];

    void setSampleRate(float sR)
    {
        sampleRate_ = sR;
        freqScale_  = 2.f / sR;
    }
    void setPole(float R, float freq);
    void setTd(float Td);
    void setT60(float T60);
    void update(float R, float freq, float Td, float T60);
    void process(float **__restrict in, float **__restrict out, int num);
};
