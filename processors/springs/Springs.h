#pragma once

#include "dsp/AllPass.h"
#include "dsp/LFO.h"
#include "dsp/LinAlg.h"
#include "dsp/MultiRate.h"
#include "dsp/Noise.h"
#include "dsp/Smoother.h"
#include "dsp/VAFilters.h"

#ifdef SPRINGS_RMS
#include "dsp/RMS.h"
#include "dsp/Stack.h"
#endif

namespace processors
{

class Springs
{
  public:
    static constexpr auto N = 4;

    static constexpr auto MaxBlockSize = 512;
    static constexpr auto MaxDecimate  = 8;
    static constexpr int LoopLength    = static_cast<int>(0.25 * 48000);
    static constexpr auto CascadeL     = 160;

    static constexpr auto DecimateMaxFreq = 0.78f;
    static constexpr auto DCBlockFreq     = 10.f;

    static constexpr auto EQPeak      = 100.f;
    static constexpr auto EQBandWidth = 5.f;

    static constexpr auto LowPassRes = 0.6f;

    static constexpr auto MinRWithMaxCascadeL = 0.1f;

    static constexpr auto NonLinearityGain = 0.2f;

    static constexpr auto APDiffMax = 0.5f;
    static constexpr auto APDiffMin = 0.08f;

    static constexpr float freqFactor[]   = {0.98f, 1.02f, 0.97f, 1.03f};
    static constexpr float RFactor[]      = {1.08f, 0.97f, 1.05f, 0.98f};
    static constexpr float loopTdFactor[] = {
        0.8293183583208989f, 1.1876863056468745f, 0.94273342f, 1.2432815625f};
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
        /* loop modulation */
        dsp::iNoise<N> noise;
        loopMod_.setPhase(noise.process());
    }

    // set processor samplerate
    template <class ReAlloc = decltype(std::realloc)>
    void prepare(float sampleRate, int blockSize,
                 ReAlloc realloc = std::realloc);

    template <class Free = decltype(std::free)>
    void free(Free free = std::free);

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
    void setFreq(float freq, int blockSize);
    void setRes(float R, int blockSize);
    void setTd(float Td, int blockSize);
    void setChaos(float chaos, int blockSize);
    void setT60(float T60, int blockSize);
    void setDiffusion(float dif, int blockSize);
    void setScatter(float scatter, int blockSize);
    void setDryWet(float drywet, int blockSize);
    void setWidth(float width, int blockSize);

    // update
    void update(float R, float freq, float Td, float T60, float diffusion,
                float chaos, float scatter, float width, float drywet,
                int blockSize);

    // main process
    void process(const float *const *__restrict in,
                 float *const *__restrict out, int num);

  private:
    int M_{1};
    float sampleRate_{48000.f};
    float freqScale_{2.f / 48000.f};
    int maxBlockSize_{};

    float drywet_{0.f};
    float width_{1.f};
    float R_{0.f};
    float freq_{0.f};
    float Td_{0.f};
    float T60_{0.f};
    float diffusion_{0.f};
    float chaos_{0.f};
    float scatter_{1.f};

    dsp::ControlSmoother<1> dry_{{1.f}};
    dsp::ControlSmoother<2> wet_{{}};

    static constexpr auto diffScatterFactor = 0.22f;
    static constexpr auto minScatter        = 0.1f;
    auto getScatterFactor() const
    {
        return minScatter + scatter_ +
               (diffusion_ * diffScatterFactor) * (1 - scatter_);
    }

    // allpass
    dsp::AllPass2<NAP> allpass_{};
    dsp::fSample<NAP> allpassIntermediary_{};
    typename decltype(allpass_)::State allpassState_[APCascadeL]{};
    unsigned int apNStages_{APCascadeL};

    dsp::va::SVF<N, dsp::va::LowPass> lowpass_{};
    typename decltype(lowpass_)::State lowpassState_{};

    dsp::va::SVF<N, dsp::va::BandPass> eq_{};
    typename decltype(eq_)::State eqState_{};

    dsp::va::OnePole<N, dsp::va::HighPass> dcblocker_{};
    typename decltype(dcblocker_)::State dcblockerState_{};

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

    dsp::DelayLine<LoopLength / 2> predelaydl_;
    dsp::TapNoInterp<N> predelay_;

    float loopGain_{0.f};
    static constexpr auto defaultTd = LoopLength * 0.1f;
    dsp::ControlSmoother<N> loopTd_{
        {defaultTd, defaultTd, defaultTd, defaultTd}};
    dsp::fData<N> loopModAmp_{};
    dsp::LFOParabolic<N> loopMod_{};
    dsp::Noise<N> loopChaosNoise_{};
    dsp::SmootherLin<N> loopChaos_{};
    dsp::fData<N> loopChaosMod_{};
    using LoopType = dsp::TapCubic<N>;
    dsp::DelayLine<LoopLength, nextTo(predelaydl_)> loopdl_;

    dsp::CopyDelayLine<N, 1, nextTo(loopdl_)> loopRippleDL_;

    dsp::AllPass<N, dsp::TapNoInterp<N>> ap1_{{}};
    dsp::DelayLine<LoopLength / 5, nextTo(loopRippleDL_)> ap1dl_;

    dsp::linalg::fMatrix<N> mixMatrix_;

    static constexpr auto BufDecSize = nextTo(ap1dl_) + MaxBlockSize;
    dsp::Buffer<dsp::fSample<N>, BufDecSize> bufferDec_;

    dsp::fSample<N> *x_{};
    dsp::fSample<N> *xdecimate_{};

// section for rms output of springs
#ifdef SPRINGS_RMS
  public:
    static constexpr auto RMSSize      = 128;
    static constexpr auto RMSOverlap   = RMSSize - 32;
    static constexpr auto RMSStackSize = 64;

    const dsp::fSample<N> *getRMSStack() const
    {
        return rmsStack_.getSamples();
    }
    const size_t *getRMSStackPos() const { return rmsStack_.getPos(); }

  private:
    dsp::RMS<N, RMSSize, RMSOverlap> rms_;
    dsp::Stack<N, RMSStackSize> rmsStack_;
#endif
};

template <class ReAlloc>
void Springs::prepare(float sampleRate, int blockSize, ReAlloc realloc)
{
    sampleRate_   = sampleRate;
    freqScale_    = 2.f / sampleRate;
    maxBlockSize_ = std::min(blockSize, MaxBlockSize);

    /* loop size modulation */
    dsp::fData<N> freq;
    for (size_t i = 0; i < N; ++i) {
        freq[i] = loopModFreq[i] * freqScale_;
    }
    loopMod_.setFreq(freq);

    // alloc ressources
    x_         = (dsp::fSample<N> *)realloc(x_, sizeof(dsp::fSample<N>) *
                                                    static_cast<size_t>(maxBlockSize_));
    xdecimate_ = (dsp::fSample<N> *)realloc(
        xdecimate_,
        sizeof(dsp::fSample<N>) * static_cast<size_t>((maxBlockSize_ + 1) / 2));

    // set buffers
#define allocateBuffer(buffer)                                \
    {                                                         \
        auto *b = buffer.getBuffer();                         \
        constexpr auto size =                                 \
            sizeof(dsp::fSample<N>) * decltype(buffer)::Size; \
        b = (dsp::fSample<N> *)realloc(b, size);              \
        memset(b, 0, size);                                   \
        buffer.setBuffer(b);                                  \
    }

    allocateBuffer(buffer_);
    allocateBuffer(bufferDec_);

#undef allocateBuffer
}

template <class Free> void Springs::free(Free free)
{
    free(x_);
    x_ = nullptr;

    free(buffer_.getBuffer());
    buffer_.setBuffer(nullptr);

    free(bufferDec_.getBuffer());
    bufferDec_.setBuffer(nullptr);
}
} // namespace processors
