#pragma once

#include "dsp/AllPass.h"
#include "dsp/LFO.h"
#include "dsp/MultiRate.h"
#include "dsp/Noise.h"
#include "dsp/Smoother.h"
#include "dsp/VAFilters.h"

namespace processors
{

class Springs
{
  public:
    static constexpr auto N = 4;

    static constexpr auto MaxBlockSize = 512;
    static constexpr auto MaxDecimate  = 8;
    static constexpr int LoopLength    = 0.25 * 48000;
    static constexpr auto CascadeL     = 160;

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

    Springs()
    {
        buffer_.setBuffer(bufferArray_);
        bufferDec_.setBuffer(bufferDecArray_);

        /* loop modulation */
        dsp::iNoise<N> noise;
        loopMod_.setPhase(noise.process());
    }

    // set processor samplerate
    void setSampleRate(float sR);

    // getters
    float getDryWet() const { return drywet_; }
    float getWidth() const { return width_; }
    float getFreq() const { return freq_; }
    float getR() const { return R_; }
    float getTd() const { return Td_; }
    float getChaos() const { return chaos_; }
    float getT60() const { return T60_; }
    float getDiffusion() const { return diffusion_; }

    // setters
    void setFreq(float R, float freq);
    void setTd(float Td, float chaos);
    void setT60(float T60);
    void setDiffusion(float dif);
    void setSpread(float spread);
    void setDryWet(float drywet) { drywet_ = drywet; }
    void setWidth(float width);

    // update
    void update(float R, float freq, float Td, float T60, float diffusion,
                float chaos, float spread, float width, float drywet);

    // main process
    void process(const float *const *__restrict in,
                 float *const *__restrict out, int num);

  private:
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
    float spread_{1.f};

    float widthcos_{1.f}, widthsin_{0.f};

    // allpass
    dsp::va::SVF<NAP, dsp::va::AllPass> allpass_;
    dsp::fSample<NAP> allpassIntermediary_{0.f};
    typename decltype(allpass_)::State allpassState_[APCascadeL];
    unsigned int apNStages_{APCascadeL};

    dsp::va::SVF<N, dsp::va::LowPass> lowpass_;
    typename decltype(lowpass_)::State lowpassState_;

    dsp::va::SVF<N, dsp::va::BandPass> eq_;
    typename decltype(eq_)::State eqState_;

    dsp::va::OnePole<N, dsp::va::HighPass> dcblocker_;
    typename decltype(dcblocker_)::State dcblockerState_;

    /* decimate and interpolate memory lines */
    using MR = dsp::MultiRate<N, 15, MaxDecimate>;
    MR::DLDecimate dldecimate_;
    MR::DLInterpolate dlinterpolate_;
    int decimateId_{0};

    static constexpr auto BufSize = nextTo(dldecimate_) + MaxBlockSize;

    // multirate pointer
  public:
    using MRs = MR::WithBuffer<BufSize>;

  private:
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
};
} // namespace processors
