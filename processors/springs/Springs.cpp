#include "Springs.h"
#include "dsp/Hadamard.h"

namespace processors
{
// multirate converter
static const Springs::MRs multirates{Springs::DecimateMaxFreq};

void Springs::setSampleRate(float sR)
{
    sampleRate_ = sR;
    freqScale_  = 2.f / sR;

    dsp::fData<N> freq;
    for (size_t i = 0; i < N; ++i) {
        freq[i] = loopModFreq[i] * freqScale_;
    }
    loopMod_.setFreq(freq);
}

void Springs::update(float R, float freq, float Td, float T60, float diffusion,
                     float chaos, float scatter, float width, float drywet)
{
    if (R != R_ || freq != freq_) {
        setFreq(R, freq);
    }
    if (Td != Td_ || chaos != chaos_) {
        setTd(Td, chaos);
    }
    if (T60 != T60_) {
        setT60(T60);
    }
    if (diffusion != diffusion_) {
        setDiffusion(diffusion);
    }
    if (width != width_) {
        setWidth(width);
    }
    if (scatter_ != scatter) {
        setScatter(scatter);
    }
    setDryWet(drywet);
}

void Springs::setFreq(float R, float freq)
{
    auto freqScaled = freq * freqScale_;
    int M           = DecimateMaxFreq / freqScaled;
    M               = std::min(M, MaxDecimate);

    freqScaled *= M;
    dsp::fData<N> freqs;
    dsp::fData<NAP> freqsAP;
    dsp::fData<NAP> Rs;
    for (size_t i = 0; i < N; ++i) {
        auto fFactor = 1.f + (freqFactor[i] - 1.f) * getScatterFactor();
        auto rFactor = 1.f + (RFactor[i] - 1.f) * getScatterFactor();
        freqs[i]     = freqScaled * fFactor;
        freqs[i]     = std::min(0.995f, freqs[i]);
        freqs[i]     = std::max(0.005f, freqs[i]);
        freqsAP[i]   = freqs[i];
        Rs[i]        = std::abs(R) * rFactor;
    }

    if (R < 0) {
        for (size_t i = 0; i < N; ++i) {
            freqsAP[i] = 1.f - freqsAP[i];
        }
    }
    for (size_t i = N; i < NAP; ++i) {
        freqsAP[i] = freqsAP[i % N];
        Rs[i]      = Rs[i % N];
    }

    allpass_.setFreq(freqsAP, Rs);

    multirate_  = multirates.get(M);
    decimateId_ = 0;

    lowpass_.setFreq(freqs);

    int oldM = M_;
    M_       = M;
    R_       = R;
    freq_    = freq;

    if (M != oldM) {

        auto eqPeak = EQPeak * freqScale_ * M;
        eq_.setFreq({eqPeak, eqPeak, eqPeak, eqPeak});
        eq_.setBandWidth({EQBandWidth, EQBandWidth, EQBandWidth, EQBandWidth});

        auto dcblockfreq = DCBlockFreq * freqScale_ * M;
        dcblocker_.setFreq(
            {dcblockfreq, dcblockfreq, dcblockfreq, dcblockfreq});

        setTd(Td_, chaos_);
    }
}

void Springs::setTd(float Td, float chaos)
{
    Td_    = Td;
    chaos_ = chaos;
    dsp::iData<N> loopEchoT;
    dsp::iData<N> predelayT;
    dsp::fSample<N> loopTd;
    float sampleTd = Td * sampleRate_ / M_;
    for (size_t i = 0; i < N; ++i) {
        auto loopFactor = 1.f + (loopTdFactor[i] - 1.f) * getScatterFactor();
        loopTd[i]       = sampleTd * loopFactor;

        loopModAmp_[i]   = loopTd[i] * loopModFactor[i];
        loopChaosMod_[i] = loopTd[i] * 0.07f * std::pow(chaos, 2.5f);

        loopEchoT[i] = loopTd[i] / 5.f;

        predelayT[i] = loopTd[i] * 0.5f;

        loopTd[i] -= loopEchoT[i];
    }

    loopTd_.set(loopTd, invBlockSize_ * M_);

    predelay_.setDelay(predelayT);
    ap1_.setDelay(loopEchoT);

    setT60(T60_);
}

void Springs::setDiffusion(float dif)
{
    diffusion_ = dif;

    auto apdif = APDiffMin + dif * (APDiffMax - APDiffMin);
    ap1_.setCoeff({-apdif, -apdif, -apdif, -apdif});

    mixMatrix_ = dsp::hadamardInterpolMatrix<N>(dif);
    setTd(Td_, chaos_);
    setFreq(R_, freq_);
}

void Springs::setT60(float T60)
{
    T60_      = T60;
    loopGain_ = -powf(0.001, Td_ / T60_);
}

void Springs::setWidth(float width)
{
    width_     = width;
    auto theta = M_PIf / 4.f * (1.f - width_);
    widthcos_  = cosf(theta);
    widthsin_  = sinf(theta);
}

void Springs::setScatter(float scatter)
{
    scatter_ = scatter;
    setTd(Td_, chaos_);
    setFreq(R_, freq_);
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

    int blockSize = std::min(count, MaxBlockSize);

    auto ctxt    = dsp::BufferContext(x_, blockSize, buffer_);
    auto ctxtdec = dsp::BufferContext(xdecimate_, blockSize, bufferDec_);

    auto noise = loopChaosNoise_.process();
    for (size_t i = 0; i < N; ++i) {
        noise[i] *= loopChaosMod_[i];
    }
    loopChaos_.set(noise, sampleRate_ * 0.05f);

    while (count) {

        contextFor(ctxt)
        {
            auto &x = c.getSignal();
            arrayFor(x, k)
            {
                auto mix = (*(inl++) + *(inr++)) * 0.5f;
                arrayFor(x[0], i) { x[k][i] = mix; }
            }
        }

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

            // mixing matrix (hadamard);
            arrayFor(loop, k) { loop[k] = mixMatrix_.mult(loop[k]); }

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
        contextFor(ctxtdec.vec())
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

            float wet[2];
            wet[0] = widthcos_ * mix[0] + widthsin_ * mix[1];
            wet[1] = widthcos_ * mix[1] + widthsin_ * mix[0];

            drywet_.step();
            auto drywet = drywet_.get()[0][0];
            *outl       = *inl2 + drywet * (wet[0] - *inl2);
            *outr       = *inr2 + drywet * (wet[1] - *inr2);
            ++outl, ++inl2;
            ++outr, ++inr2;
        }

        ctxt.nextBlock();
        ctxtdec.nextBlock();

        count -= blockSize;
    }
    ctxt.save(buffer_);
    ctxtdec.save(bufferDec_);

    drywet_.reset();
    loopTd_.reset();
}
} // namespace processors
