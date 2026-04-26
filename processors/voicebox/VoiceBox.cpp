#include "VoiceBox.h"
#include <iostream>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

static inline float atune(float freq)
{
    auto m  = 12.f * std::log2(freq / (440.f / 24000.f));
    auto im = std::round(m / 3) * 3;
    return 440.f / 24000.f * std::pow(2.f, (im) / 12.f);
}

static inline auto barkScale(float freq)
{
    return -1.0674f * std::sqrt(2.f / dsp::constants<float>::pi *
                                std::atan(0.06583f * freq / 1000.f)) +
           0.1916f;
}

void VoiceBox::prepare(float sampleRate)
{
    pitchTracker_.prepare(sampleRate);

    auto freqScale = 2.f / sampleRate;

    sampleRate_ = sampleRate;
    zffAnalysis_.setWindowSize(
        static_cast<int>(sampleRate * kAnalysisWindowSize));

    filter_.setWarpIn(barkScale(sampleRate));
    filter_.setWarpOut(barkScale(sampleRate));

    /* signal.cheby2(8,80,1,"high", analog=True, output="sos") */
    static constexpr decltype(dcblocker_)::SOS kSosAnalog = {
        {{{1.f, 0.f, 0.03806023f, 1.f, 2.15704486f, 1.24729461f},
          {1.f, 0.f, 0.30865828f, 1.f, 1.82865434f, 1.51789266f},
          {1.f, 0.f, 0.69134172f, 1.f, 1.82865434f, 1.51789266f},
          {1.f, 0.f, 0.96193977f, 1.f, 1.82865434f, 1.51789266f}}}};
    dcblocker_.fromAnalog(kSosAnalog, 40.f * freqScale);

    grainGenerator_.setFreq(100.f * freqScale);
    grainGenerator_.start();
}

void VoiceBox::update(float shift, float formantShift)
{
    if (shift != getPitchShift()) {
        setPitchShift(shift);
    }
    if (formantShift != getFormantShift()) {
        setFormantShift(formantShift);
    }
}

void VoiceBox::setPitchShift(float shift)
{
    shift_ = shift;
    psola_.setShift(shift);
}

void VoiceBox::setFormantShift(float formantShift)
{
    formantShift_ = formantShift;
    psola_.setPlayRate(formantShift);
    grainGenerator_.setPlayRate(formantShift);
    filter_.setWarpOut(barkScale(sampleRate_ / formantShift));
}

void VoiceBox::setRange(PitchTracker::Range range)
{
    pitchTracker_.setRange(range);
}

void VoiceBox::process(const float *const *__restrict in,
                       float *const *__restrict out, int count)
{
    const float *min = *in;
    float *mout      = *out;

    if (min != mout) {
        for (int i = 0; i < count; ++i) mout[i] = min[i];
    }

    auto ctxt = dsp::Context(mout, count);
    CTXTRUN(ctxt) { dcblocker_.process(ctxt, dcblockerState_); };

    auto pitch  = pitchBuffer_;
    auto voiced = voicedBuffer_;
    while (count) {
        auto blockSize = std::min(kMaxBlockSize, count);

        float *pitchTrackOut[] = {pitch, voiced};
        pitchTracker_.process(in, pitchTrackOut, blockSize);

        // set pitch atomic value to last pitch computed
        auto storeFreq = pitch[blockSize - 1] * sampleRate_ * 0.5f;
        storeFreq      = voiced[blockSize - 1] > 0.95f ? storeFreq : 0.f;
        atomicFreq_.store(storeFreq);

        auto ctxt       = dsp::BufferContext(mout, blockSize, buffer_);
        auto ctxtPitch  = dsp::Context(pitch, blockSize);
        auto ctxtVoiced = dsp::Context(voiced, blockSize);

        CTXTRUN(ctxtPitch) { auto p = ctxtPitch.getInput(); };

        int n = 0;
        CTXTRUN2(ctxt, ctxtPitch)
        {
            zffAnalysis_.process(ctxt, analysisDL_, [&, this](auto &marker) {
                // ctxt.setOutput(ctxt.getInput()-1.f);
                auto freq = ctxtPitch.getInput();
                assert(marker.period_ >= 0.f);

                if (voiced[n] > 0.9f) {
                    // ctxt.setOutput(ctxt.getInput()+2.f);
                    auto trackedPeriod = 2.f / freq;
                    while (marker.period_ > 1.1f * trackedPeriod) {
                        auto ratio  = marker.period_ / trackedPeriod;
                        auto iratio = std::round(ratio);
                        float newPeriod;
                        auto lastMarker = zffAnalysis_.getMarker();
                        if (std::abs(iratio - ratio) > 0.2) {
                            newPeriod = lastMarker.period_;
                        } else {
                            newPeriod = marker.period_ / iratio;
                            newPeriod += 0.8f * (trackedPeriod - newPeriod);
                        }
                        std::cout << "oops \t" << trackedPeriod << " \t"
                                  << marker.period_ << " \t" << ratio << " \t"
                                  << iratio << " \t" << newPeriod << "\n";
                        assert(newPeriod > 0.f);
                        auto newPosition =
                            ctxt.modAdd(lastMarker.position_, newPeriod);
                        zffAnalysis_.addMarker({newPosition, newPeriod});
                        marker.period_ -= trackedPeriod;
                    }
                    auto lastMarker = zffAnalysis_.getMarker();
                    if (marker.period_ > 0.9f * trackedPeriod) {
                        marker.period_ +=
                            0.8f * (trackedPeriod - marker.period_);
                        marker.position_ =
                            ctxt.modAdd(lastMarker.position_, marker.period_);
                        assert(marker.period_ > 0.f);

                        // ctxt.setOutput(ctxt.getInput()+2.f);
                        return true;
                    }
                    return false;
                }
                return true;
            });
            ++n;
        };

        n = 0;
        CTXTRUN2(ctxt, ctxtPitch)
        {
            static float v = 0.f;
            voiced[n]      = v += 0.002f * (voiced[n] - v);
            if (ctxtPitch.getInput() > 1e-6f) {
                psola_.setFreq(ctxtPitch.getInput());
                // psola_.setFreq(atune(ctxtPitch.getInput() * shift_));
                // psola_.setFreq(400.f/24000.f*shift_);
            }
            if (!voiced_ && voiced[n] > 0.2f) {
                // set period of last marker
                auto period = 2.f / ctxtPitch.getInput();
                zffAnalysis_.getClosestMarker(ctxt).period_ = period;

                // start psola synthesis
                psola_.start(ctxt, zffAnalysis_);
                voiced_ = true;
            } else if (voiced_ && voiced[n] < 0.2f) {
                // stop psola synthesis
                psola_.stop(ctxt);
                // zffAnalysis_.resetMarkers();
                voiced_ = false;
            }
            if (voiced[n] > 0.99f) {
                grainGenerator_.stop();
            }
            if (voiced[n] < 0.6f) {
                grainGenerator_.start();
            }
            auto x = ctxt.getInput();

            // x += noise_.process()*0.002f;
            // rls_.process(dsp::Context{&x}, warpAState1_, filter_);
            // auto cpyFilter = filter_;
            // cpyFilter.smoothCoeffs(0.92f);

            // cpyFilter.analyze(ctxt, warpAState_);
            // cpyFilter.compensateResidualAnalyze(ctxt, warpAState_);

            // cpyFilter.compensateResidualReconstruct(ctxt, warpRState_);
            // cpyFilter.reconstruct(ctxt, warpRState_);

            psola_.process(ctxt, synthDL_, zffAnalysis_);
            grainGenerator_.process(dsp::BufferContext(&x, 1, ctxt.getBuffer()),
                                    synthDL_);
            x += (ctxt.getInput() - x) * voiced[n];
            ctxt.setOutput(x);

            ++n;
        };

        count -= blockSize;
        buffer_.nextBlock(ctxt);
        mout += blockSize;
        pitch += blockSize;
        voiced += blockSize;
    }
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
