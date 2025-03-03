

#include "../../AllPass.h"
#include "../../LFO.h"
#include "../../MultiRate.h"
#include "../../Noise.h"
#include "../../Smoother.h"
#include "../../VAFilters.h"

static constexpr auto N = 4;

class Springs
{
  public:
    static constexpr auto MaxBlockSize = 512;
    static constexpr auto MaxDecimate  = 8;
    static constexpr int LoopLength    = 0.25 * 48000;
    static constexpr auto CascadeL     = 120;

    static constexpr auto DecimateMaxFreq = 0.78f;
    static constexpr auto DCBlockFreq     = 10.f;

    static constexpr auto EQPeak      = 100.f;
    static constexpr auto EQBandWidth = 5.f;

    static constexpr auto NonLinearityGain = 0.2f;

    static constexpr float freqFactor[]    = {0.98f, 1.02f, 0.97f, 1.03f};
    static constexpr float RFactor[]       = {1.03f, 0.97f, 1.05f, 0.98f};
    static constexpr float loopTdFactor[]  = {0.979f, 1.0f, 1.035f, 1.05f};
    static constexpr float loopModFreq[]   = {0.2f, 0.4f, 0.2f, 0.3f};
    static constexpr float loopModFactor[] = {0.0045f, 0.003f, 0.005f, 0.0037f};
    static constexpr float loopRippleGain  = 0.016f;

    // use simd to improve all pass chain
    // for example
    //  [x1,x2,x3,x4] -----> AP[4]*L -----> [y1,y2,y3,y4]
    // becomes
    //  [x1,x2,x3,x4,y1,y2,y3,y4] --> AP[8]*L/2 --> [y1,y2,y3,y4,z1,z2,z3,z4]
    // with y as intermediary values and z as final values
    // reducing the number of computation by 2
    // introduce a delay in chain of size APChainSize
    static constexpr auto APChainSize = SIMDSIZE / sizeof(float) / N;
    static constexpr auto NAP         = N * APChainSize;
    static constexpr auto APCascadeL  = CascadeL / APChainSize;

    Springs(float sampleRate)
    {
        buffer_.setBuffer(bufferArray_);
        bufferDec_.setBuffer(bufferDecArray_);
        setSampleRate(sampleRate);

        /* loop modulation */
        dsp::iNoise<N> noise;
        loopMod_.setPhase(noise.process());

        dsp::fData<N> freq;
        for (int i = 0; i < N; ++i) {
            freq[i] = loopModFreq[i] * freqScale_;
        }
        loopMod_.setFreq(freq);
    }

    int M_{1};
    float sampleRate_{48000.f};
    float freqScale_{2.f / 48000.f};

    float drywet_{1.f};
    float width_{1.f};
    float R_{0.f};
    float freq_{0.f};
    float Td_{0.f};
    float T60_{0.f};
    float diffusion_{0.f};
    float chaos_{0.f};

    float widthcos_{1.f}, widthsin_{0.f};

    dsp::va::SVF<NAP, dsp::va::AllPass> allpass_;
    dsp::fSample<NAP> allpassIntermediary_{0.f};
    typename decltype(allpass_)::State allpassState_[APCascadeL]{{{{0.f}}}};

    dsp::va::SVF<N, dsp::va::LowPass> lowpass_;
    typename decltype(lowpass_)::State lowpassState_;

    dsp::va::SVF<N, dsp::va::BandPass> eq_;
    typename decltype(eq_)::State eqState_{{{0.f}}};

    dsp::va::OnePole<N, dsp::va::HighPass> dcblocker_;
    typename decltype(dcblocker_)::State dcblockerState_{{0.f}};

    using MR = dsp::MultiRate<N, 15, MaxDecimate>;
    MR::DLDecimate dldecimate_;
    MR::DLInterpolate dlinterpolate_;
    int decimateId_{0};

    static constexpr auto BufSize = nextTo(dldecimate_) + MaxBlockSize;

    using MRs = MR::WithBuffer<BufSize>;
    const MRs::Base *multirate_;

    dsp::Buffer<dsp::fSample<N>, BufSize> buffer_;
    dsp::fSample<N> bufferArray_[decltype(buffer_)::Size]{{0.f}};

    dsp::DelayLine<LoopLength / 2> predelaydl_;
    dsp::TapNoInterp<N> predelay_;

    float loopGain_{0.f};
    dsp::fData<N> loopTd_{0.f};
    dsp::fData<N> loopModAmp_{0.f};
    dsp::LFOParabolic<N> loopMod_;
    dsp::Noise<N> loopChaosNoise_;
    dsp::SmootherLin<N> loopChaos_;
    dsp::fData<N> loopChaosMod_;
    using LoopType = dsp::TapCubic<N>;
    dsp::DelayLine<LoopLength, nextTo(predelaydl_)> loopdl_;

    dsp::CopyDelayLine<N, 1, nextTo(loopdl_)> loopRippleDL_;

    dsp::AllPass<N, dsp::TapNoInterp<N>> ap1_{{0.f, 0.f, 0.f, 0.f}};
    dsp::DelayLine<LoopLength / 5, nextTo(loopRippleDL_)> ap1dl_;

    static constexpr auto BufDecSize = nextTo(ap1dl_) + MaxBlockSize;
    dsp::Buffer<dsp::fSample<N>, BufDecSize> bufferDec_;
    dsp::fSample<N> bufferDecArray_[decltype(bufferDec_)::Size]{{0.f}};

    dsp::fSample<N> x_[MaxBlockSize]{{0.f}};
    dsp::fSample<N> xdecimate_[MaxBlockSize / 2]{{0.f}};

    void setSampleRate(float sR)
    {
        sampleRate_ = sR;
        freqScale_  = 2.f / sR;
    }

    void setFreq(float R, float freq);
    void setTd(float Td, float chaos);
    void setT60(float T60);
    void setDiffusion(float dif);
    void update(float R, float freq, float Td, float T60, float diffusion,
                float chaos, float width, float drywet);
    void process(float **__restrict in, float **__restrict out, int num);
};
