#include "Springs.h"
#include "dsp/Orthogonal.h"

namespace processors
{
// multirate converter
static const Springs::MRs multirates{Springs::DecimateMaxFreq};

void Springs::update(float R, float freq, float Td, float T60, float tone,
                     float chaos, float scatter, float width, float drywet,
                     int blockSize)
{
    if (freq != freq_) {
        setFreq(freq, blockSize);
    }
    if (R != R_) {
        setRes(R, blockSize);
    }
    if (Td != Td_) {
        setTd(Td, blockSize);
    }
    if (chaos != chaos_) {
        setChaos(chaos, blockSize);
    }
    if (T60 != T60_) {
        setT60(T60, blockSize);
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
    auto M          = static_cast<int>(DecimateMaxFreq / freqScaled);
    M               = std::min(M, MaxDecimate);

    freqScaled *= static_cast<float>(M);
    dsp::fData<N> freqs;
    dsp::fData<NAP> freqsAP;
    for (size_t i = 0; i < N; ++i) {
        auto fFactor = 1.f + (freqFactor[i] - 1.f) * getScatterFactor();
        freqs[i]     = freqScaled * fFactor;
        freqs[i]     = std::min(0.995f, freqs[i]);
        freqs[i]     = std::max(0.005f, freqs[i]);
        freqsAP[i]   = freqs[i];
    }

    if (R_ < 0) {
        for (size_t i = 0; i < N; ++i) {
            freqsAP[i] = 1.f - freqsAP[i];
        }
    }
    for (size_t i = N; i < NAP; ++i) {
        freqsAP[i] = freqsAP[i % N];
    }

    allpass_.setFreq(freqsAP);
    lowpass_.setFreq(freqs, {
                                LowPassRes,
                                LowPassRes,
                                LowPassRes,
                                LowPassRes,
                            });

    multirate_  = multirates.get(M);
    decimateId_ = 0;

    int oldM = M_;
    M_       = M;
    freq_    = freq;

    if (M != oldM) {
        auto fM = static_cast<float>(M);

        auto dcblockfreq = DCBlockFreq * freqScale_ * fM;
        dcblocker_.setFreq(
            {dcblockfreq, dcblockfreq, dcblockfreq, dcblockfreq});

        setTd(Td_, blockSize);
        setTone(tone_, blockSize);
    }
}

void Springs::setRes(float R, int /*blockSize*/)
{
    R_ = R;

    dsp::fData<NAP> Rs;
    for (size_t i = 0; i < N; ++i) {
        auto rFactor = 1.f + (RFactor[i] - 1.f) * getScatterFactor();
        Rs[i]        = std::abs(R_) * rFactor;
    }

    for (size_t i = N; i < NAP; ++i) {
        Rs[i] = Rs[i % N];
    }

    allpass_.setRes(Rs);

    /* if abs(R) smaller than a certain value, reduce the cascade size
     * this helps to avoid long ringing around allpass phasing frequency */
    if (std::abs(R) < MinRWithMaxCascadeL) {
        apNStages_ = static_cast<unsigned int>(
            std::abs(R) / MinRWithMaxCascadeL * APCascadeL);
    } else {
        apNStages_ = APCascadeL;
    }
}

void Springs::setTd(float Td, int blockSize)
{
    Td_ = Td;
    dsp::iData<N> loopEchoT;
    dsp::iData<N> predelayT;
    dsp::fSample<N> loopTd;
    float sampleTd = Td * sampleRate_ / static_cast<float>(M_);
    for (size_t i = 0; i < N; ++i) {
        auto loopFactor = 1.f + (loopTdFactor[i] - 1.f) * getScatterFactor();
        loopTd[i]       = sampleTd * loopFactor;

        loopModAmp_[i]   = loopTd[i] * loopModFactor[i];
        loopChaosMod_[i] = loopTd[i] * 0.07f * std::pow(chaos_, 2.5f);

        loopEchoT[i] = static_cast<int>(loopTd[i] / 5.f);

        predelayT[i] = static_cast<int>(loopTd[i] * .5f);

        loopTd[i] -= static_cast<float>(loopEchoT[i]);
    }

    loopTd_.set(loopTd, static_cast<float>(M_) / static_cast<float>(blockSize));

    predelay_.setDelay(predelayT);
    ap1_.setDelay(loopEchoT);

    setT60(T60_, blockSize);
}

void Springs::setChaos(float chaos, int blockSize)
{
    chaos_ = chaos;
    setTd(Td_, blockSize);
}

void Springs::setTone(float tone, int /*blockSize*/)
{
    tone_ = tone;

    auto eqPeak = dsp::expScale(ToneMin, ToneMax, tone) * freqScale_ * M_;
    static constexpr auto maxEqPeak = 0.95f;
    eqPeak                          = std::min(maxEqPeak, eqPeak);
    eq_.setFreq({eqPeak, eqPeak, eqPeak, eqPeak});
    eq_.setBandWidth({EQBandWidth, EQBandWidth, EQBandWidth, EQBandWidth});
}

void Springs::setT60(float T60, int /*blockSize*/)
{
    T60_      = T60;
    loopGain_ = -powf(0.001f, Td_ / T60_);
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
    setTd(Td_, blockSize);
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

            // allpass diffusion
            {
                auto apctxt = c;
                apctxt.setSamples(loop);
                ap1_.process(apctxt, ap1dl_);
            }

            // householder transform
            arrayFor(loop, k) { loop[k] = dsp::householder(loop[k]); }

            loopRippleDL_.write(c, loop);
            auto loopRipple = dsp::TapTail{}.read(c, loopRippleDL_);

            inFor(loop, k, i)
            {
                loopRipple[k][i] +=
                    loopRippleGain * (loop[k][i] - loopRipple[k][i]);
            }

            inFor(x, k, i) { x[k][i] += loopRipple[k][i] * loopGain_; }
        }

        contextFor(ctxtdec)
        {
            auto &x = c.getSignal();
            inFor(x, k, i)
            {
                x[k][i] *= NonLinearityGain;
                x[k][i] = dsp::tanh(x[k][i]);
                x[k][i] /= NonLinearityGain;
            }
        }
        contextFor(ctxtdec) { dcblocker_.process(c, dcblockerState_); }

        // allpass cascade
        auto allpassIntermediary = allpassIntermediary_;
        contextFor(ctxtdec)
        {
            auto &x = c.getSignal();

            // shift intermediary values
            for (size_t j = APChainSize - 1; j > 0; --j) {
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

        contextFor(ctxtdec) { lowpass_.process(c, lowpassState_); }
        contextFor(ctxtdec)
        {
            auto &x = c.getSignal();
            loopdl_.write(c, x);
        }

        contextFor(ctxtdec) { eq_.process(c, eqState_); }

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
