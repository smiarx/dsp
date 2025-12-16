#include "VoiceBox.h"

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

void VoiceBox::prepare(float sampleRate)
{
    constexpr auto kAmpDecay = 0.015f;
    ampFollower_.setRate(0.1f, static_cast<float>(kDecFactor) /
                                   (kAmpDecay * sampleRate));
}

void VoiceBox::update(float forgetFactor, float warpIn, float warpOut,
                      float shift)
{
    rls_.setForgetFactor(forgetFactor);
    filter_.setWarpIn(warpIn);
    filter_.setWarpOut(warpOut);
    pitchShift_.setShift(shift);
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

    CTXTRUN(ctxtDec)
    {
        auto amp = dsp::abs(ctxtDec.getInput());
        ampFollower_.process(dsp::Context(&amp), ampFollowerState_);

        auto origFilter = filter_;
        if (amp > 0.01f) {
            pitchTracker_.process(ctxtDec, pitchTrackerState_);
            auto freq = pitchTrackerState_.getFreq();
            pitchShift_.setFreq(freq);
            rls_.process(ctxtDec, warpAState_, filter_);
        }
        filter_.compensateResidualAnalyze(ctxtDec, warpAState_);

        pitchShift_.process(ctxtDec, pitchShiftDL_);

        filter_.compensateResidualReconstruct(ctxtDec, warpRState_);
        origFilter.reconstruct(ctxtDec, warpRState_);
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
