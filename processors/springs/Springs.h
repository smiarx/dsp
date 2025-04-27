#pragma once

#include "dsp/AllPass.h"
#include "dsp/IIRFilter.h"
#include "dsp/LFO.h"
#include "dsp/MultiRate.h"
#include "dsp/Noise.h"
#include "dsp/Smoother.h"
#include "dsp/VAFilters.h"

#ifdef SPRINGS_RMS
#include "dsp/RMS.h"
#include "dsp/Stack.h"
#endif

#ifdef SPRINGS_SHAKE
#include "dsp/Enveloppe.h"
#endif

namespace processors
{

class Springs
{
  public:
    static constexpr auto N = 4;

    static constexpr auto kDefaultSamplerRate = 48000.f;

    static constexpr auto kMaxBlockSize         = 512;
    static constexpr auto kMaxDecimate          = 8;
    static constexpr auto kMaxLoopLengthSeconds = 0.2f;
    static constexpr auto kLoopLengthSeconds    = 0.25f;
    static constexpr auto kLoopLength =
        static_cast<int>(kLoopLengthSeconds * kDefaultSamplerRate);
    static constexpr auto kCascadeL = 160;

    static constexpr auto kDecimateMaxFreq = 0.78f;
    static constexpr auto kDcBlockFreq     = 10.f;

    static constexpr auto kEqBandWidth = 5.f;

    static constexpr auto kMinRWithMaxCascadeL = 0.1f;

    static constexpr auto kNonLinearityGain = 0.2f;

    static constexpr auto kToneMin = 80.f;
    static constexpr auto kToneMax = 5000.f;

    static constexpr float kFreqFactor[]   = {0.98f, 1.02f, 0.97f, 1.03f};
    static constexpr float kRFactor[]      = {1.08f, 0.97f, 1.05f, 0.98f};
    static constexpr float kLoopTdFactor[] = {
        0.8293183583208989f, 1.1876863056468745f, 0.94273342f, 1.2432815625f};
    static constexpr float kLoopModFreq[]   = {0.2f, 0.4f, 0.2f, 0.3f};
    static constexpr float kLoopModFactor[] = {0.0035f, 0.002f, 0.0038f,
                                               0.0027f};
    static constexpr float kLoopRippleGain  = 0.016f;

    // use simd to improve all pass chain
    // for example
    //  [x1,x2,x3,x4] -----> AP[4]*L -----> [y1,y2,y3,y4]
    // becomes
    //  [x1,x2,x3,x4,y1,y2,y3,y4] --> AP[8]*L/2 --> [y1,y2,y3,y4,z1,z2,z3,z4]
    // with y as intermediary values and z as final values
    // reducing the number of computation by 2
    // introduce a delay in chain of size APChainSize
    static constexpr auto kApChainSize = SIMDSIZE / sizeof(float) / N;
    static constexpr auto kNap         = N * kApChainSize;
    static constexpr auto kApCascadeL  = kCascadeL / kApChainSize;

    Springs()
    {
        /* loop modulation */
        dsp::INoise<N> noise;
        loopMod_.setPhase(noise.process());
    }

    // set processor samplerate
    template <class ReAlloc = decltype(std::realloc)>
    void prepare(float sampleRate, int blockSize,
                 ReAlloc realloc = std::realloc);

    template <class Free = decltype(std::free)>
    void free(Free free = std::free);

    // getters
    [[nodiscard]] float getDryWet() const { return drywet_; }
    [[nodiscard]] float getWidth() const { return width_; }
    [[nodiscard]] float getFreq() const { return freq_; }
    [[nodiscard]] float getR() const { return r_; }
    [[nodiscard]] float getTd() const { return td_; }
    [[nodiscard]] float getChaos() const { return chaos_; }
    [[nodiscard]] float getT60() const { return t60_; }
    [[nodiscard]] float getTone() const { return tone_; }

    // setters
    void setFreq(float freq, int blockSize);
    void setRes(float r, int blockSize);
    void setTd(float td, int blockSize);
    void setChaos(float chaos, int blockSize);
    void setT60(float t60, int blockSize);
    void setTone(float tone, int blockSize);
    void setScatter(float scatter, int blockSize);
    void setDryWet(float drywet, int blockSize);
    void setWidth(float width, int blockSize);

    // update
    void update(float r, float freq, float td, float t60, float tone,
                float chaos, float scatter, float width, float drywet,
                int blockSize);

    // main process
    void process(const float *const *__restrict in,
                 float *const *__restrict out, int count);

  private:
    int rateFactor_{1};
    float sampleRate_{kDefaultSamplerRate};
    float freqScale_{2.f / kDefaultSamplerRate};
    int maxBlockSize_{};

    float drywet_{0.f};
    float width_{1.f};
    float r_{0.f};
    float freq_{0.f};
    float td_{0.f};
    float t60_{0.f};
    float tone_{0.f};
    float chaos_{0.f};
    float scatter_{1.f};

    dsp::ControlSmoother<1> dry_{{1.f}};
    dsp::ControlSmoother<2> wet_{{}};

    static constexpr auto kMinScatter = 0.1f;
    [[nodiscard]] auto getScatterFactor() const
    {
        return kMinScatter + scatter_;
    }

    // allpass
    dsp::AllPass2<kNap> allpass_{};
    dsp::fSample<kNap> allpassIntermediary_{};
    typename decltype(allpass_)::State allpassState_[kApCascadeL]{};
    unsigned int apNStages_{kApCascadeL};

    dsp::IIRFilter<N, 10> lowpass_{};
    typename decltype(lowpass_)::Mem<N> lowpassState_{};

    dsp::va::SVF<N, dsp::va::kBandPass> eq_{};
    typename decltype(eq_)::State eqState_{};

    dsp::va::OnePole<N, dsp::va::kHighPass> dcblocker_{};
    typename decltype(dcblocker_)::State dcblockerState_{};

    /* decimate and interpolate memory lines */
    using MR = dsp::MultiRate<N, 15, kMaxDecimate>;
    MR::DLDecimate dldecimate_;
    MR::DLInterpolate dlinterpolate_;
    int decimateId_{0};

    static constexpr auto kBufSize = nextTo(dldecimate_) + kMaxBlockSize;

    // multirate pointer
  public:
    using MRs = MR::WithBuffer<kBufSize>;

  private:
    const MRs::Base *multirate_;

    dsp::Buffer<dsp::fSample<N>, kBufSize> buffer_;

    dsp::DelayLine<kLoopLength / 2> predelaydl_;
    dsp::TapNoInterp<N> predelay_;

    float loopGain_{0.f};
    static constexpr auto kDefaultTd = kLoopLength * 0.1f;
    dsp::ControlSmoother<N> loopTd_{
        {kDefaultTd, kDefaultTd, kDefaultTd, kDefaultTd}};
    dsp::fData<N> loopModAmp_{};
    dsp::LFOParabolic<N> loopMod_{};
    dsp::Noise<N> loopChaosNoise_{};
    dsp::SmootherLin<N> loopChaos_{};
    dsp::fData<N> loopChaosMod_{};
    using LoopType = dsp::TapCubic<N>;
    dsp::DelayLine<kLoopLength, nextTo(predelaydl_)> loopdl_;

    dsp::CopyDelayLine<N, 1, nextTo(loopdl_)> loopRippleDL_;

    static constexpr auto kBufDecSize = nextTo(loopRippleDL_) + kMaxBlockSize;
    dsp::Buffer<dsp::fSample<N>, kBufDecSize> bufferDec_;

    dsp::fSample<N> *x_{};
    dsp::fSample<N> *xdecimate_{};

// section for rms output of springs
#ifdef SPRINGS_RMS
  public:
    static constexpr auto kRmsSize      = 128;
    static constexpr auto kRmsOverlap   = kRmsSize - 32;
    static constexpr auto kRmsStackSize = 64;

    [[nodiscard]] const dsp::fSample<N> *getRMSStack() const
    {
        return rmsStack_.getSamples();
    }
    [[nodiscard]] const size_t *getRMSStackPos() const
    {
        return rmsStack_.getPos();
    }

  private:
    dsp::RMS<N, kRmsSize, kRmsOverlap> rms_;
    dsp::Stack<N, kRmsStackSize> rmsStack_;
#endif

#ifdef SPRINGS_SHAKE
  public:
    static constexpr auto kShakeEnvUp   = 0.001f;
    static constexpr auto kShakeEnvDown = 0.045f;
    void shake()
    {
        shakeEnv_.set({1.f}, {kShakeEnvUp * sampleRate_},
                      {kShakeEnvDown * sampleRate_});
    }

  private:
    dsp::DoubleRamp<1> shakeEnv_;
    dsp::Noise<1> shakeNoise_;
#endif
};

template <class ReAlloc>
void Springs::prepare(float sampleRate, int blockSize, ReAlloc realloc)
{
    sampleRate_   = sampleRate;
    freqScale_    = 2.f / sampleRate;
    maxBlockSize_ = std::min(blockSize, kMaxBlockSize);

    /* loop size modulation */
    dsp::fData<N> freq;
    for (size_t i = 0; i < N; ++i) {
        freq[i] = kLoopModFreq[i] * freqScale_;
    }
    loopMod_.setFreq(freq);

    // alloc ressources
    x_         = (dsp::fSample<N> *)realloc(x_, sizeof(dsp::fSample<N>) *
                                                    static_cast<size_t>(maxBlockSize_));
    xdecimate_ = (dsp::fSample<N> *)realloc(
        xdecimate_,
        sizeof(dsp::fSample<N>) * static_cast<size_t>((maxBlockSize_ + 1) / 2));

    // set buffers
#define allocateBuffer(buffer)                                 \
    {                                                          \
        auto *b = buffer.getBuffer();                          \
        constexpr auto size =                                  \
            sizeof(dsp::fSample<N>) * decltype(buffer)::kSize; \
        b = (dsp::fSample<N> *)realloc(b, size);               \
        memset(b, 0, size);                                    \
        buffer.setBuffer(b);                                   \
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
