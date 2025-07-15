#include "Springs.h"
#include "dsp/Buffer.h"
#include "dsp/Context.h"
#include "dsp/FastMath.h"
#include "dsp/Orthogonal.h"
#include "dsp/Utils.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>

namespace processors
{
// multirate converter
const Springs::MRD Springs::kDecimate{};
const Springs::MRI Springs::kInterpolate{};

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

    auto fFactor = 1.f + (kFreqFactor.load() - 1.f) * getScatterFactor();
    auto freqs   = freqScaled * fFactor;
    freqs        = dsp::min(mtype(0.995f).load(), freqs);
    freqs        = dsp::max(mtype(0.005f).load(), freqs);
    auto freqsAP = dsp::batch<mtype>::simdtype::convert(freqs);

    if (r_ < 0) {
        freqsAP /= 2;
    }

    allpass_.setFreq(freqsAP);

    /* scipy.signal.cheby1(10,2,1,analog=True,output='sos') */
    decltype(lowpass_)::tf kTfAnalog = {{
        {{{0., 0., 0.00255383}, {1., 0.21436212, 0.0362477}}},
        {{{0., 0., 1.}, {1., 0.19337886, 0.21788333}}},
        {{{0., 0., 1.}, {1., 0.15346633, 0.51177596}}},
        {{{0., 0., 1.}, {1., 0.09853145, 0.80566858}}},
        {{{0., 0., 1.}, {1., 0.03395162, 0.98730422}}},
    }};
    lowpass_.tfAnalog(kTfAnalog, freqs);

    decimateId_ = 0;

    int oldRateFactor = rateFactor_;
    rateFactor_       = rateFactor;
    freq_             = freq;

    if (rateFactor != oldRateFactor) {
        auto fM = static_cast<float>(rateFactor);

        auto dcblockfreq = kDcBlockFreq * freqScale_ * fM;
        dcblocker_.setFreq(dcblockfreq);

        setTd(td_, blockSize);
        setTone(tone_, blockSize);
    }
}

void Springs::setRes(float r, int blockSize)
{
    bool signChanged = std::signbit(r) != std::signbit(r_);
    r_               = r;

    auto rFactor = 1.f + (kRFactor.load() - 1.f) * getScatterFactor();
    auto rs      = dsp::abs(r_) * rFactor;

    allpass_.setRes(dsp::batch<mtype>::simdtype::convert(rs));

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
    td_            = td;
    float sampleTd = td * sampleRate_ / static_cast<float>(rateFactor_);

    // prevent from going higher than buffer size
    sampleTd = std::min(sampleTd, kMaxLoopLengthSeconds * kDefaultSamplerRate *
                                      kDefaultSamplerRate / sampleRate_);

    auto loopFactor = 1.f + (kLoopTdFactor.load() - 1.f) * getScatterFactor();
    auto loopTd     = sampleTd * loopFactor;
    loopModAmp_     = loopTd * kLoopModFactor;
    loopChaosMod_   = loopTd * 0.07f * dsp::pow(chaos_, 2.5f);
    auto predelayT  = dsp::toInt(loopTd * .5f);

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
    eq_.setFreq(eqPeak);
    eq_.setBandWidth(kEqBandWidth);
}

void Springs::setT60(float t60, int /*blockSize*/)
{
    t60_      = t60;
    loopGain_ = dsp::pow(0.001f, td_ / t60_);
}

void Springs::setWidth(float width, int blockSize)
{
    width_          = width;
    auto theta      = dsp::constants<float>::pi_4 * (1.f - width_);
    auto wet        = dsp::sin(dsp::constants<float>::pi_2 * drywet_);
    auto wetchannel = dsp::cos(theta) * wet;
    auto wetcross   = dsp::sin(theta) * wet;

    float invBlockSize = 1.f / static_cast<float>(blockSize);
    wet_.set({wetchannel, wetcross},
             invBlockSize * static_cast<float>(rateFactor_));
}

void Springs::setDryWet(float drywet, int blockSize)
{
    drywet_            = drywet;
    float invBlockSize = 1.f / static_cast<float>(blockSize);
    float dry          = dsp::cos(dsp::constants<float>::pi_2 * drywet);
    dry_.set(dry, invBlockSize);

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
    noise *= loopChaosMod_;
    loopChaos_.set(noise, static_cast<int>(sampleRate_ * 0.05f));

    while (count) {

        blockSize = std::min(count, blockSize);

        auto ctxtIn =
            dsp::BufferContext(reinterpret_cast<float *>(x_) +
                                   (static_cast<ptrdiff_t>(blockSize * 3)),
                               blockSize, bufferIn_);
        auto ctxtOut =
            dsp::Context(reinterpret_cast<dsp::mfloat<2> *>(x_), blockSize);
        auto ctxtdec = dsp::BufferContext(x_, blockSize, bufferDec_);

        CTXTRUN(ctxtIn)
        {
            auto mix = (*(inl++) + *(inr++)) * 0.5f;
            ctxtIn.setOutput(mix);
        };

#ifdef SPRINGS_SHAKE
        // shake springs
        {
            auto ctxt = ctxtIn;
            for (int n = 0; n < ctxt.getBlockSize(); ++n, ctxt.next()) {
                auto x = ctxt.getInput();
                if (shakeEnv_.isRunning()) {
                    auto env   = shakeEnv_.process();
                    auto noise = shakeNoise_.process();
                    auto shake = env * noise;
                    x += shake;
                    ctxt.setOutput(x);
                } else {
                    break;
                }
            }
        }
#endif

        kDecimate.decimate(rateFactor_, ctxtIn, ctxtdec, dldecimate_,
                           decimateId_);

        CTXTRUN(ctxtdec)
        {
            auto x = ctxtdec.getInput();
            predelaydl_.write(ctxtdec, x);

            x = predelay_.read(ctxtdec, predelaydl_);
            ctxtdec.setOutput(x);
        };

        CTXTRUN(ctxtdec)
        {
            auto x = ctxtdec.getInput();

            LoopType looptap;
            auto mod   = loopMod_.process();
            auto chaos = loopChaos_.step();

            auto delay = loopTd_.step(ctxtdec);

            delay += mod * loopModAmp_ + chaos;
            looptap.setDelay(delay);

            auto loop = looptap.read(ctxtdec, loopdl_);

            // householder transform
            loop = dsp::householder(loop);

            loopRippleDL_.write(ctxtdec, loop);
            auto loopRipple = dsp::TapTail{}.read(ctxtdec, loopRippleDL_);

            loopRipple += kLoopRippleGain * (loop - loopRipple);

            x += loopRipple * loopGain_;
            ctxtdec.setOutput(x);
        };

        CTXTRUN(ctxtdec)
        {
            auto x = ctxtdec.getInput();
            x *= kNonLinearityGain;
            x = dsp::tanh(x);
            x /= kNonLinearityGain;
            ctxtdec.setOutput(x);
        };
        CTXTRUN(ctxtdec) { dcblocker_.process(ctxtdec, dcblockerState_); };

        // allpass cascade
        CTXTRUN(ctxtdec)
        {
            auto x = ctxtdec.getInput();

            // shift intermediary values
            for (size_t j = kApChainSize - 1; j > 0; --j) {
                allpassIntermediary_[j] = allpassIntermediary_[j - 1];
            }
            // set value as first entries of intermediary values
            allpassIntermediary_[0] = x;

            // compute allpass filters
            dsp::Context c1(
                reinterpret_cast<dsp::batch<mtype> *>(&allpassIntermediary_));
            for (size_t j = 0; j < apNStages_; ++j) {
                allpass_.process(c1, allpassState_[j]);
            }

            // outout is last intermediary value
            x = allpassIntermediary_[kApChainSize - 1];
            ctxtdec.setOutput(x);
        };

        CTXTRUN(ctxtdec)
        {
            auto x = ctxtdec.getInput();
            loopdl_.write(ctxtdec, x);
        };

        CTXTRUN(ctxtdec) { eq_.process(ctxtdec, eqState_); };

        CTXTRUN(ctxtdec) { lowpass_.process(ctxtdec, lowpassState_); };

#ifdef SPRINGS_RMS
        rms_.processBlock(ctxtdec, rmsStack_);
#endif

        auto ctxtdecOut = dsp::BufferContext(
            reinterpret_cast<dsp::mfloat<2> *>(&x_[blockSize]) -
                (static_cast<ptrdiff_t>(ctxtdec.getBlockSize())),
            ctxtdec.getBlockSize(), bufferOut_);
        CTXTRUNREV2(ctxtdec, ctxtdecOut)
        {
            ctxtdecOut.setOutput(dsp::reduce<2>(ctxtdec.getInput()));
        };
        CTXTRUN(ctxtdecOut)
        {
            auto mix = ctxtdecOut.getInput();

            mix *= 2.f / kN;

            dsp::mfloat<2> wetsig;
            auto wet  = wet_.step(ctxtOut);
            wetsig[0] = dsp::sum(wet * mix);
            wetsig[1] = dsp::sum(wet * dsp::flip<1>(mix));

            ctxtdecOut.setOutput(wetsig);
        };

        decimateId_ = kInterpolate.interpolate(rateFactor_, ctxtdecOut, ctxtOut,
                                               dlinterpolate_, decimateId_);

        CTXTRUN(ctxtOut)
        {
            auto wetsig = ctxtOut.getInput();
            auto dry    = dry_.step(ctxtOut);
            *outl       = *inl2 * dry + wetsig[0];
            *outr       = *inr2 * dry + wetsig[1];
            ++outl, ++inl2;
            ++outr, ++inr2;
        };

        bufferIn_.nextBlock(ctxtIn);
        bufferOut_.nextBlock(ctxtdecOut);
        bufferDec_.nextBlock(ctxtdec);

        count -= blockSize;
    }

    dry_.reset();
    wet_.reset();
    loopTd_.reset();
}
} // namespace processors
