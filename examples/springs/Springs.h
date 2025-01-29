

#include "../../AllPass.h"
#include "../../MultiRate.h"
#include "../../Filter.h"

namespace _dsp = dsp;

static constexpr auto N = 4;

class Springs
{
public:
    Springs (float sampleRate) {
        buffer_.setBuffer(bufferArray_);
        bufferDec_.setBuffer(bufferDecArray_);
        setSampleRate(sampleRate);
    }

    static constexpr auto MaxBlockSize = 512;
    static constexpr auto MaxDecimate = 8;
    static constexpr auto CascadeL = 120;

    int M_{1};
    float sampleRate_{48000.f};
    float freqScale_{1.f/48000.f};
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

    static constexpr auto BufSize = nextTo(dldecimate_)+MaxBlockSize;

    using MRs = MR::WithBuffer<BufSize>;
    const MRs::Base* multirate_;

    dsp::Buffer<dsp::Signal<N>, BufSize> buffer_;
    dsp::Signal<N> bufferArray_[decltype(buffer_)::Size] = {{0.f}};

    static constexpr int LoopLength = 0.25*48000;
    float loopGain_{0.f};
    dsp::TapAllPass<N, dsp::TapNoInterp<N>> looptap_;
    //dsp::TapNoInterp<N> looptap_;
    dsp::NestedDelayLine<dsp::DelayLine<LoopLength>, dsp::CopyDelayLine<N,1>> loopdl_;
    //dsp::DelayLine<LoopLength> loopdl_;

    static constexpr auto BufDecSize = nextTo(loopdl_)+MaxBlockSize;
    dsp::Buffer<dsp::Signal<N>, BufDecSize> bufferDec_;
    dsp::Signal<N> bufferDecArray_[decltype(bufferDec_)::Size] = {{0.f}};

    dsp::Signal<N> x_[MaxBlockSize];
    dsp::Signal<N> xdecimate_[MaxBlockSize/2];

    void setSampleRate(float sR)
    {
        sampleRate_ = sR;
        freqScale_ = 2.f/sR;
    }
    void setPole(float R, float freq);
    void setTd(float Td);
    void setT60(float T60);
    void update(float R, float freq, float Td, float T60);
    void process(float**__restrict in, float**__restrict out, int num);
};
