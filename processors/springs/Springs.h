#pragma once

#include "dsp/AllPass.h"
#include "dsp/Buffer.h"
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
    static constexpr auto kN = 4;
    using type               = float;
    using mtype              = dsp::multi<type, kN>;

    static constexpr type kDefaultSamplerRate = 48000.;

    static constexpr auto kMaxBlockSize         = 512;
    static constexpr auto kMaxDecimate          = 8;
    static constexpr type kMaxLoopLengthSeconds = 0.2;
    static constexpr type kLoopLengthSeconds    = 0.25;
    static constexpr auto kLoopLength =
        static_cast<int>(kLoopLengthSeconds * kDefaultSamplerRate);
    static constexpr auto kCascadeL = 160;

    static constexpr type kDecimateMaxFreq = 0.78;
    static constexpr type kDcBlockFreq     = 10.;

    static constexpr type kEqBandWidth = 5.;

    static constexpr type kMinRWithMaxCascadeL = 0.1;

    static constexpr type kNonLinearityGain = 0.2;

    static constexpr type kToneMin = 80.;
    static constexpr type kToneMax = 5000.;

    static constexpr mtype kFreqFactor   = {0.98, 1.02, 0.97, 1.03};
    static constexpr mtype kRFactor      = {1.08, 0.97, 1.05, 0.98};
    static constexpr mtype kLoopTdFactor = {
        0.8293183583208989, 1.1876863056468745, 0.94273342, 1.2432815625};
    static constexpr mtype kLoopModFreq   = {0.2, 0.4, 0.2, 0.3f};
    static constexpr mtype kLoopModFactor = {0.0035, 0.002, 0.0038, 0.0027};
    static constexpr type kLoopRippleGain = 0.016;

    // use simd to improve all pass chain
    // for example
    //  [x1,x2,x3,x4] -----> AP[4]*L -----> [y1,y2,y3,y4]
    // becomes
    //  [x1,x2,x3,x4,y1,y2,y3,y4] --> AP[8]*L/2 --> [y1,y2,y3,y4,z1,z2,z3,z4]
    // with y as intermediary values and z as final values
    // reducing the number of computation by 2
    // introduce a delay in chain of size APChainSize
    static constexpr auto kApChainSize = DSP_MAX_VEC_SIZE / sizeof(type) / kN;
    static constexpr auto kApCascadeL  = kCascadeL / kApChainSize;

    Springs()
    {
        /* loop modulation */
        dsp::INoise<mtype> noise;
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
    type sampleRate_{kDefaultSamplerRate};
    type freqScale_{2.f / kDefaultSamplerRate};
    int maxBlockSize_{};

    type drywet_{0.f};
    type width_{1.f};
    type r_{0.f};
    type freq_{0.f};
    type td_{0.f};
    type t60_{0.f};
    type tone_{0.f};
    type chaos_{0.f};
    type scatter_{1.f};

    dsp::ControlSmoother<type> dry_{1.f};
    dsp::ControlSmoother<dsp::mfloat<2>> wet_{};

    static constexpr type kMinScatter = 0.1;
    [[nodiscard]] auto getScatterFactor() const
    {
        return kMinScatter + scatter_;
    }

    // allpass
    dsp::AllPass2<dsp::batch<mtype>> allpass_{};
    std::array<mtype, kApChainSize> allpassIntermediary_{};
    typename decltype(allpass_)::State allpassState_[kApCascadeL]{};
    unsigned int apNStages_{kApCascadeL};

    dsp::IIRFilter<mtype, 10> lowpass_{};
    typename decltype(lowpass_)::State lowpassState_{};

    dsp::va::SVF<mtype, dsp::va::kBandPass> eq_{};
    typename decltype(eq_)::State eqState_{};

    dsp::va::OnePole<mtype, dsp::va::kHighPass> dcblocker_{};
    typename decltype(dcblocker_)::State dcblockerState_{};

    /* decimate and interpolate memory lines */
    using MRD = dsp::MultiRateDecimate<float, 15, kMaxDecimate>;
    using MRI = dsp::MultiRateInterpolate<mtype, 15, kMaxDecimate>;
    static const MRD kDecimate;
    static const MRI kInterpolate;
    MRD::DLDecimate<0> dldecimate_;
    MRI::DLInterpolate<0> dlinterpolate_;
    int decimateId_{0};

    static constexpr auto kBufSize = nextTo(dldecimate_) + kMaxBlockSize;
    dsp::Buffer<type, kBufSize> buffer_;

    dsp::DelayLine<kLoopLength / 2, nextTo(dlinterpolate_)> predelaydl_;
    dsp::TapNoInterp<mtype> predelay_;

    type loopGain_{};
    static constexpr type kDefaultTd    = kLoopLength * 0.1;
    dsp::ControlSmoother<mtype> loopTd_ = {
        {kDefaultTd, kDefaultTd, kDefaultTd, kDefaultTd}};
    mtype loopModAmp_{};
    dsp::lfo::Parabolic<mtype> loopMod_{};
    dsp::Noise<mtype> loopChaosNoise_{};
    dsp::SmootherLin<mtype> loopChaos_{};
    mtype loopChaosMod_{};
    using LoopType = dsp::TapCubic<mtype>;
    dsp::DelayLine<kLoopLength, nextTo(predelaydl_)> loopdl_;

    dsp::CopyDelayLine<mtype, 1, nextTo(loopdl_)> loopRippleDL_;

    static constexpr auto kBufDecSize = nextTo(loopRippleDL_) + kMaxBlockSize;
    dsp::Buffer<mtype, kBufDecSize> bufferDec_;

    mtype *x_{};
    mtype *xdecimate_{};

// section for rms output of springs
#ifdef SPRINGS_RMS
  public:
    static constexpr auto kRmsSize      = 128;
    static constexpr auto kRmsOverlap   = kRmsSize - 32;
    static constexpr auto kRmsStackSize = 64;

    [[nodiscard]] const auto *getRMSStack() const
    {
        return rmsStack_.getSamples();
    }
    [[nodiscard]] const size_t *getRMSStackPos() const
    {
        return rmsStack_.getPos();
    }

  private:
    dsp::RMS<mtype, kRmsSize, kRmsOverlap> rms_;
    dsp::Stack<mtype, kRmsStackSize> rmsStack_;
#endif

#ifdef SPRINGS_SHAKE
  public:
    static constexpr auto kShakeEnvUp   = 0.001f;
    static constexpr auto kShakeEnvDown = 0.045f;
    void shake()
    {
        shakeEnv_.set(1.f, {kShakeEnvUp * sampleRate_},
                      {kShakeEnvDown * sampleRate_});
    }

  private:
    dsp::DoubleRamp<type> shakeEnv_;
    dsp::Noise<type> shakeNoise_;
#endif
};

template <class ReAlloc>
void Springs::prepare(float sampleRate, int blockSize, ReAlloc realloc)
{
    sampleRate_   = sampleRate;
    freqScale_    = 2.f / sampleRate;
    maxBlockSize_ = std::min(blockSize, kMaxBlockSize);

    /* loop size modulation */
    mtype freq = dsp::load(kLoopModFreq) * freqScale_;
    loopMod_.setFreq(freq);

    // alloc ressources
    x_         = (mtype *)realloc(x_,
                                  sizeof(mtype) * static_cast<size_t>(maxBlockSize_));
    xdecimate_ = (mtype *)realloc(
        xdecimate_,
        sizeof(mtype) * static_cast<size_t>((maxBlockSize_ + 1) / 2));

    // set buffers
#define allocateBuffer(buffer, type)                                   \
    {                                                                  \
        auto *b             = buffer.getData();                        \
        constexpr auto size = sizeof(mtype) * decltype(buffer)::kSize; \
        b                   = (type *)realloc(b, size);                \
        memset(b, 0, size);                                            \
        buffer.setData(b);                                             \
    }

    allocateBuffer(buffer_, float);
    allocateBuffer(bufferDec_, mtype);

#undef allocateBuffer
}

template <class Free> void Springs::free(Free free)
{
    free(x_);
    x_ = nullptr;
    free(xdecimate_);
    xdecimate_ = nullptr;

    free(buffer_.getData());
    buffer_.setData(nullptr);

    free(bufferDec_.getData());
    bufferDec_.setData(nullptr);
}
} // namespace processors
