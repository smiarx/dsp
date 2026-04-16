#include "VoiceBox.h"
#include <climits>
#include <iostream>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

static inline float atune(float freq)
{
    auto m  = 12.f * std::log2(freq / (440.f / 24000.f));
    auto im = std::round(m / 6) * 6;
    return std::pow(2.f, (im - m) / 12.f);
}

void VoiceBox::prepare(float sampleRate)
{
    constexpr auto kAmpDecay = 0.002f;
    ampFollower_.setRate(0.1f, static_cast<float>(kDecFactor) /
                                   (kAmpDecay * sampleRate));

    // pitchShift_.setFreq(0.01f);

    decltype(ptLowpass_)::SOS kTfAnalog = {{{
        {0.f, 0.f, 0.01021531f, 1.f, 0.26637237f, 0.05650064f},
        {0.f, 0.f, 1.f, 1.f, 0.22581959f, 0.32709869f},
        {0.f, 0.f, 1.f, 1.f, 0.15088783f, 0.70978212f},
        {0.f, 0.f, 1.f, 1.f, 0.05298476f, 0.98038017f},
    }}};
    ptLowpass_.fromAnalog(kTfAnalog, kDecFactor * 2 * 600 / sampleRate);

    zcOp_.setCutoff(2 * 10 / sampleRate);
}

void VoiceBox::update(float forgetFactor, float warpIn, float warpOut,
                      float ptop, float ptrho, float ptlambda, float shift)
{
    rls_.setForgetFactor(forgetFactor);
    filter_.setWarpIn(warpIn);
    filter_.setWarpOut(warpOut);
    ptOp1_.setCoeff(ptop);
    ptOp2_.setCoeff(ptop);
    fastPitchTracker_.setRho(ptrho);
    fastPitchTracker_.setLambda(ptlambda);
    psola_.setShift(shift);
    shift_ = shift;
}

void VoiceBox::process(const float *const *__restrict in,
                       float *const *__restrict out, int count)
{
    const float *min = *in;
    float *mout      = *out;

    if (min != mout) {
        for (int i = 0; i < count; ++i) mout[i] = min[i];
    }

    auto ctxt    = dsp::BufferContext(mout, count, buffer_);
    auto ctxtDec = dsp::BufferContext(mout, count, bufferDec_);

    decimate_.decimate(ctxt, ctxtDec, decimateDL_, decimateId_);
    auto decBS = ctxtDec.getBlockSize();

    // dsp::Noise<float> noise;
    CTXTRUN(ctxtDec)
    {

        // auto origFilter = filter_;
        //  auto zcValue = zc_.process(ctxtDec) ? 1.f : 0.f;
        //  zcOp_.process(dsp::Context(&zcValue), zcOpState_);

        auto x = ctxtDec.getInput();
        ptLowpass_.process(ctxtDec, ptLowpassState_);
        ptOp1_.process(ctxtDec, ptOp1St_);

        auto amp = dsp::abs(ctxtDec.getInput());
        amp *= amp;
        ampFollower_.process(dsp::Context(&amp), ampFollowerState_);
        // ptOp2_.process(ctxtDec, ptOp2St_);
        // if (zcValue < 0.2) {
        // pitchTracker_.process(ctxtDec, pitchTrackerState_);
        // ctxtDec.setOutput(ctxtDec.getInput()*1000.f);
        pitchTracker_.process(ctxtDec, pitchTrackerDL_);

        // auto freq = pitchTrackerState_.getFreq();
        // if (amp > 0.005f) {

        float period = 0.f;
        if (chooseFastPitchTracker_) {
            fastPitchTracker_.process(ctxtDec, fastPitchTrackerState_);
            // std::cout << "fast:" << freq*12000 << "\n";
            //

            auto freq = fastPitchTrackerState_.getFreq();
            freq      = dsp::max(freq, 50.f / 12000.f);
            period    = 2.f / freq;
            assert(period > 0);
            // ctxtDec.setOutput(freq);
            //
            // auto shift = atune(freq);
            // std::cout << "freq " << freq*12000.f << " - shift " << shift <<
            // "\n";
            psola_.setShift(shift_);

            auto pole = std::abs(fastPitchTrackerState_.getPole());
            pole *= pole;
            ampFollower_.process(dsp::Context(&pole), poleFollowerState_);

            if (poleFollowerState_ < ampFollowerState_ * 10000.f) {
                std::cout << "fast stop\n";
                chooseFastPitchTracker_ = false;
                psola_.stop(ctxtDec);
            } else {
                //    //std::cout << "fast:" << freq*24000 << "\n";
                //    //pitchShift_.setFreq(freq);
                //    //ctxtDec.setOutput(freq);
            }
        }

        if (!chooseFastPitchTracker_) {
            period = pitchTracker_.findPeak(0.95f);
            assert(period >= 0);
            if (period > 0.f) {

                auto freq = 2.f / period;
                freq      = dsp::max(freq, 50.f / 24000.f);

                std::cout << "slow:" << freq * 12000 << "\n";
                // std::cout << "slow:" << period << "\n";
                chooseFastPitchTracker_ = true;
                poleFollowerState_      = 5000;
                fastPitchTrackerState_.setFreq(freq);

                psola_.start(ctxtDec, period, pitchShiftDL_);
            } else {
                // pitchShift_.setShift(1);
                //  ctxtDec.setOutput(0);
            }
        }

        // pitchShift_.setShift(shift_);
        // ctxtDec.setOutput(freq);
        //}
        // else
        //{
        //    pitchShift_.setShift(1);
        //}

        ctxtDec.setOutput(x);

        rls_.process(ctxtDec, warpAState_, filter_);

        //} else {
        //    //filter_.analyze(ctxtDec, warpAState_);
        //}

        filter_.compensateResidualAnalyze(ctxtDec, warpAState_);

        psola_.findMarkers(ctxtDec, pitchShiftDL_, period);
        psola_.process(ctxtDec, pitchShiftDL_);

        // ctxtDec.setOutput(noise.process()*0.1f);

        filter_.compensateResidualReconstruct(ctxtDec, warpRState_);
        // origFilter.reconstruct(ctxtDec, warpRState_);
    };

    // move data to end of buffer to avoid write conflict
    for (int i = 0; i < decBS; ++i) {
        mout[count - i - 1] = mout[decBS - i - 1];
    }
    auto ctxtInt = dsp::BufferContext(&mout[count - decBS], decBS, bufferDec_);
    decimateId_ =
        interpolate_.interpolate(ctxtInt, ctxt, interpolateDL_, decimateId_);

    buffer_.nextBlock(ctxt);
    bufferDec_.nextBlock(ctxtDec);
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
