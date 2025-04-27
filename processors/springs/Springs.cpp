#include "Springs.h"
#include "dsp/Buffer.h"
#include "dsp/Context.h"
#include "dsp/FastMath.h"
#include "dsp/Orthogonal.h"
#include "dsp/Signal.h"
#include "dsp/Utils.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>

namespace processors
{
// multirate converter
static const Springs::MRs kMultirates{Springs::kDecimateMaxFreq};

void Springs::update(float r, float freq, float td, float t60, float tone,
                     float chaos, float scatter, float width, float drywet,
                     int blockSize)
{
    if (freq != freq_) {
        setFreq(freq, blockSize);
    }
    if (r != r_) {
        setRes(r, blockSize);
    }
    if (td != td_) {
        setTd(td, blockSize);
    }
    if (chaos != chaos_) {
        setChaos(chaos, blockSize);
    }
    if (t60 != t60_) {
        setT60(t60, blockSize);
    }
    if (tone != tone_) {
        setTone(tone, blockSize);
    }
    if (width != width_) {
        setWidth(width, blockSize);
    }
    if (scatter_ != scatter) {
        setScatter(scatter, blockSize);
    }
    if (drywet != drywet_) {
        setDryWet(drywet, blockSize);
    }
}

void Springs::setFreq(float freq, int blockSize)
{
    auto freqScaled = freq * freqScale_;
    auto rateFactor = static_cast<int>(kDecimateMaxFreq / freqScaled);
    rateFactor      = std::min(rateFactor, kMaxDecimate);

    freqScaled *= static_cast<float>(rateFactor);
    dsp::fData<N> freqs;
    dsp::fData<kNap> freqsAP;
    for (size_t i = 0; i < N; ++i) {
        auto fFactor = 1.f + (kFreqFactor[i] - 1.f) * getScatterFactor();
        freqs[i]     = freqScaled * fFactor;
        freqs[i]     = std::min(0.995f, freqs[i]);
        freqs[i]     = std::max(0.005f, freqs[i]);
        freqsAP[i]   = freqs[i];
    }

    if (r_ < 0) {
        for (size_t i = 0; i < N; ++i) {
            freqsAP[i] = freqsAP[i] / 2.f;
        }
    }
    for (size_t i = N; i < kNap; ++i) {
        freqsAP[i] = freqsAP[i % N];
    }

    allpass_.setFreq(freqsAP);

    /* scipy.signal.cheby1(10,2,1,analog=True,output='sos') */
    constexpr float kTfAnalog[][2][3] = {
        {{0., 0., 0.00255383}, {1., 0.21436212, 0.0362477}},
        {{0., 0., 1.}, {1., 0.19337886, 0.21788333}},
        {{0., 0., 1.}, {1., 0.15346633, 0.51177596}},
        {{0., 0., 1.}, {1., 0.09853145, 0.80566858}},
        {{0., 0., 1.}, {1., 0.03395162, 0.98730422}},
    };
    lowpass_.tfAnalog(kTfAnalog, freqs);

    multirate_  = kMultirates.get(rateFactor);
    decimateId_ = 0;

    int oldRateFactor = rateFactor_;
    rateFactor_       = rateFactor;
    freq_             = freq;

    if (rateFactor != oldRateFactor) {
        auto fM = static_cast<float>(rateFactor);

        auto dcblockfreq = kDcBlockFreq * freqScale_ * fM;
        dcblocker_.setFreq(
            {dcblockfreq, dcblockfreq, dcblockfreq, dcblockfreq});

        setTd(td_, blockSize);
        setTone(tone_, blockSize);
    }
}

void Springs::setRes(float r, int blockSize)
{
    bool signChanged = std::signbit(r) != std::signbit(r_);
    r_               = r;

    dsp::fData<kNap> rs;
    for (size_t i = 0; i < N; ++i) {
        auto rFactor = 1.f + (kRFactor[i] - 1.f) * getScatterFactor();
        rs[i]        = std::abs(r_) * rFactor;
    }

    for (size_t i = N; i < kNap; ++i) {
        rs[i] = rs[i % N];
    }

    allpass_.setRes(rs);

    /* if abs(R) smaller than a certain value, reduce the cascade size
     * this helps to avoid long ringing around allpass phasing frequency */
    if (std::abs(r) < kMinRWithMaxCascadeL) {
        apNStages_ = static_cast<unsigned int>(
            std::abs(r) / kMinRWithMaxCascadeL * kApCascadeL);
    } else {
        apNStages_ = kApCascadeL;
    }

    if (signChanged) {
        setFreq(freq_, blockSize);
    }
}

void Springs::setTd(float td, int blockSize)
{
    td_ = td;
    dsp::iData<N> predelayT;
    dsp::fSample<N> loopTd;
    float sampleTd = td * sampleRate_ / static_cast<float>(rateFactor_);

    // prevent from going higher than buffer size
    sampleTd = std::min(sampleTd, kMaxLoopLengthSeconds * kDefaultSamplerRate *
                                      kDefaultSamplerRate / sampleRate_);

    for (size_t i = 0; i < N; ++i) {
        auto loopFactor = 1.f + (kLoopTdFactor[i] - 1.f) * getScatterFactor();
        loopTd[i]       = sampleTd * loopFactor;

        loopModAmp_[i]   = loopTd[i] * kLoopModFactor[i];
        loopChaosMod_[i] = loopTd[i] * 0.07f * std::pow(chaos_, 2.5f);
        predelayT[i]     = static_cast<int>(loopTd[i] * .5f);
    }

    loopTd_.set(loopTd, static_cast<float>(rateFactor_) /
                            static_cast<float>(blockSize));

    predelay_.setDelay(predelayT);

    setT60(t60_, blockSize);
}

void Springs::setChaos(float chaos, int blockSize)
{
    chaos_ = chaos;
    setTd(td_, blockSize);
}

void Springs::setTone(float tone, int /*blockSize*/)
{
    tone_ = tone;

    auto eqPeak =
        dsp::expScale(kToneMin, kToneMax, tone) * freqScale_ * rateFactor_;
    static constexpr auto kMaxEqPeak = 0.95f;
    eqPeak                           = std::min(kMaxEqPeak, eqPeak);
    eq_.setFreq({eqPeak, eqPeak, eqPeak, eqPeak});
    eq_.setBandWidth({kEqBandWidth, kEqBandWidth, kEqBandWidth, kEqBandWidth});
}

void Springs::setT60(float t60, int /*blockSize*/)
{
    t60_      = t60;
    loopGain_ = -powf(0.001f, td_ / t60_);
}

void Springs::setWidth(float width, int blockSize)
{
    width_          = width;
    auto theta      = dsp::constants<float>::pi_4 * (1.f - width_);
    auto wet        = std::sin(dsp::constants<float>::pi_2 * drywet_);
    auto wetchannel = std::cos(theta) * wet;
    auto wetcross   = std::sin(theta) * wet;

    float invBlockSize = 1.f / static_cast<float>(blockSize);
    wet_.set({wetchannel, wetcross}, invBlockSize);
}

void Springs::setDryWet(float drywet, int blockSize)
{
    drywet_            = drywet;
    float invBlockSize = 1.f / static_cast<float>(blockSize);
    float dry          = std::cos(dsp::constants<float>::pi_2 * drywet);
    dry_.set({dry}, invBlockSize);

    setWidth(width_, blockSize);
}

void Springs::setScatter(float scatter, int blockSize)
{
    scatter_ = scatter;
    setTd(td_, blockSize);
    setFreq(freq_, blockSize);
}

void Springs::process(const float *const *__restrict in,
                      float *const *__restrict out, int count)
{
    auto *inl  = in[0];
    auto *inr  = in[1];
    auto *inl2 = inl;
    auto *inr2 = inr;
    auto *outl = out[0];
    auto *outr = out[1];

    int blockSize = maxBlockSize_;

    auto noise = loopChaosNoise_.process();
    for (size_t i = 0; i < N; ++i) {
        noise[i] *= loopChaosMod_[i];
    }
    loopChaos_.set(noise, static_cast<int>(sampleRate_ * 0.05f));

    while (count) {

        blockSize = std::min(count, blockSize);

        auto ctxt    = dsp::BufferContext(x_, blockSize, buffer_);
        auto ctxtdec = dsp::BufferContext(xdecimate_, blockSize, bufferDec_);

        contextFor(ctxt)
        {
            auto &x = c.getSignal();
            arrayFor(x, k)
            {
                auto mix = (*(inl++) + *(inr++)) * 0.5f;
                arrayFor(x[0], i) { x[k][i] = mix; }
            }
        }

#ifdef SPRINGS_SHAKE
        // shake springs
        contextFor(ctxt)
        {
            auto &x = c.getSignal();
            if (shakeEnv_.isRunning()) {
                arrayFor(x, k)
                {
                    auto env   = shakeEnv_.process();
                    auto noise = shakeNoise_.process();
                    auto shake = env[0] * noise[0];
                    arrayFor(x[0], i) { x[k][i] += shake; }
                }
            } else {
                break;
            }
        }
#endif

        multirate_->decimate(ctxt, ctxtdec, dldecimate_, decimateId_);

        contextFor(ctxtdec)
        {
            auto &x = c.getSignal();
            predelaydl_.write(c, x);

            x = predelay_.read(c, predelaydl_);

            // arrayFor(x, k) { x[k] = dsp::hadamard(x[k]); }
        }

        contextFor(ctxtdec)
        {
            auto &x = c.getSignal();

            LoopType looptap;
            auto mod   = loopMod_.process();
            auto chaos = loopChaos_.step();

            loopTd_.step();
            auto delay = loopTd_.get()[0];

            for (size_t i = 0; i < N; ++i) {
                delay[i] += mod[i] * loopModAmp_[i] + chaos[i];
            }
            looptap.setDelay(delay);

            auto loop = looptap.read(c, loopdl_);

            // householder transform
            arrayFor(loop, k) { loop[k] = dsp::householder(loop[k]); }

            loopRippleDL_.write(c, loop);
            auto loopRipple = dsp::TapTail{}.read(c, loopRippleDL_);

            inFor(loop, k, i)
            {
                loopRipple[k][i] +=
                    kLoopRippleGain * (loop[k][i] - loopRipple[k][i]);
            }

            inFor(x, k, i) { x[k][i] += loopRipple[k][i] * loopGain_; }
        }

        contextFor(ctxtdec)
        {
            auto &x = c.getSignal();
            inFor(x, k, i)
            {
                x[k][i] *= kNonLinearityGain;
                x[k][i] = dsp::tanh(x[k][i]);
                x[k][i] /= kNonLinearityGain;
            }
        }
        contextFor(ctxtdec) { dcblocker_.process(c, dcblockerState_); }

        // allpass cascade
        auto allpassIntermediary = allpassIntermediary_;
        contextFor(ctxtdec)
        {
            auto &x = c.getSignal();

            // shift intermediary values
            for (size_t j = kApChainSize - 1; j > 0; --j) {
                for (size_t i = 0; i < N; ++i) {
                    allpassIntermediary[j * N + i] =
                        allpassIntermediary[(j - 1) * N + i];
                }
            }
            // set value as first entries of intermediary values
            for (size_t i = 0; i < N; ++i) {
                allpassIntermediary[i] = x[0][i];
            }

            // compute allpass filters
            for (size_t j = 0; j < apNStages_; ++j) {
                dsp::Context c1(&allpassIntermediary);
                allpass_.process(c1, allpassState_[j]);
            }

            // outout is last intermediary value
            for (size_t i = 0; i < N; ++i) {
                x[0][i] = allpassIntermediary[i];
            }
        }
        // save intermediary values
        allpassIntermediary_ = allpassIntermediary;

        contextFor(ctxtdec)
        {
            auto &x = c.getSignal();
            loopdl_.write(c, x);
        }

        contextFor(ctxtdec) { eq_.process(c, eqState_); }

        contextFor(ctxtdec) { lowpass_.process(c, lowpassState_); }

        decimateId_ =
            multirate_->interpolate(ctxtdec, ctxt, dlinterpolate_, decimateId_);

        contextFor(ctxt)
        {
            auto &x      = c.getSignal();
            float mix[2] = {0.f};
            for (size_t i = 0; i < N / 2; ++i) {
                mix[0] += x[0][i];
            }
            for (size_t i = N / 2; i < N; ++i) {
                mix[1] += x[0][i];
            }
            mix[0] *= 2.f / N;
            mix[1] *= 2.f / N;

            float wetsig[2];
            wet_.step();
            auto wet  = wet_.get()[0];
            wetsig[0] = wet[0] * mix[0] + wet[1] * mix[1];
            wetsig[1] = wet[0] * mix[1] + wet[1] * mix[0];

            dry_.step();
            auto dry = dry_.get()[0][0];
            *outl    = *inl2 * dry + wetsig[0];
            *outr    = *inr2 * dry + wetsig[1];
            ++outl, ++inl2;
            ++outr, ++inr2;
        }

#ifdef SPRINGS_RMS
        rms_.processBlock(ctxt, rmsStack_);
#endif

        buffer_.nextBlock(ctxt);
        bufferDec_.nextBlock(ctxtdec);

        count -= blockSize;
    }

    dry_.reset();
    wet_.reset();
    loopTd_.reset();
}
} // namespace processors
