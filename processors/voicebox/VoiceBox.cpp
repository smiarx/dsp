#include "VoiceBox.h"
#include <dsp/utils/Debug.hpp>
#include <ranges>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

void VoiceBox::prepare(Hz sampleRate)
{
    sampleRate_.set(sampleRate);

    pitchTracker_.prepare(sampleRate);

    zffAnalysis_.setWindowSize(
        static_cast<int>(sampleRate_.secsToSamples(kAnalysisWindowSize).get()));

    /* signal.cheby2(8,80,1,"high", analog=True, output="sos") */
    static constexpr decltype(dcblocker_)::SOS kSosAnalog = {
        {{{1.f, 0.f, 0.03806023f, 1.f, 2.15704486f, 1.24729461f},
          {1.f, 0.f, 0.30865828f, 1.f, 1.82865434f, 1.51789266f},
          {1.f, 0.f, 0.69134172f, 1.f, 1.82865434f, 1.51789266f},
          {1.f, 0.f, 0.96193977f, 1.f, 1.82865434f, 1.51789266f}}}};
    dcblocker_.fromAnalog(kSosAnalog, sampleRate_.hzToNyquistRatio(40._Hz));

    grainGenerator_.setFreq(sampleRate_.hzToNyquistRatio(100._Hzf).get());
    grainGenerator_.start();

    voicedSmooth_.setConvergenceTime(0.01f, sampleRate_.secsToSamples(0.05_sf));

    setQuantizeSpeed(quantizeSpeed_);
    setQuantizeConfirmTime(quantizeConfirmTime_);
    setA4Frequency(a4Frequency_);
    setScale(scale_);
}

void VoiceBox::update(float shift, float formantShift, seconds quantizeSpeed,
                      seconds quantizeConfirmTime, Hz a4Frequency)
{
    if (shift != getPitchShift()) {
        setPitchShift(shift);
    }
    if (formantShift != getFormantShift()) {
        setFormantShift(formantShift);
    }
    if (quantizeSpeed != quantizeSpeed_) {
        setQuantizeSpeed(quantizeSpeed);
    }
    if (quantizeConfirmTime != quantizeConfirmTime_) {
        setQuantizeConfirmTime(quantizeConfirmTime);
    }
    if (a4Frequency != a4Frequency_) {
        setA4Frequency(a4Frequency);
    }
}

void VoiceBox::setPitchShift(float shift) { shift_ = shift; }

void VoiceBox::setFormantShift(float formantShift)
{
    formantShift_ = formantShift;
    psola_.setPlayRate(formantShift);
    grainGenerator_.setPlayRate(formantShift);
}

void VoiceBox::setQuantizeSpeed(seconds time)
{
    quantizeSpeed_ = time;
    quantizer_.setQuantizeSpeed(sampleRate_.secsToSamples(time));
}

void VoiceBox::setQuantizeConfirmTime(seconds time)
{
    quantizeConfirmTime_ = time;
    quantizer_.setQuantizeConfirmTime(sampleRate_.secsToSamples(time));
}

void VoiceBox::setA4Frequency(Hz freq)
{
    a4Frequency_      = freq;
    auto nyquistRatio = sampleRate_.hzToNyquistRatio(freq);
    pitchRef_         = dsp::PitchReference{nyquistRatio};
}

void VoiceBox::setScale(dsp::Scale scale)
{
    scale_ = scale;
    quantizer_.setScale(scale);
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
        auto lastFreq  = dsp::freq::NyquistRatio(pitch[blockSize - 1]);
        auto storeFreq = sampleRate_.nyquistRatioToHz(lastFreq);
        storeFreq      = voiced[blockSize - 1] > 0.95f ? storeFreq : 0.f;
        atomicFreq_.store(storeFreq.get());

        auto ctxt      = dsp::BufferContext(mout, blockSize, buffer_);
        auto ctxtPitch = dsp::Context(pitch, blockSize);
        // auto ctxtVoiced = dsp::Context(voiced, blockSize);

        int n = 0;
        CTXTRUN2(ctxt, ctxtPitch)
        {
            zffAnalysis_.process(ctxt, analysisDL_, [&, this](auto &marker) {
                auto freq = ctxtPitch.getInput();
                assert(marker.period_ >= 0.f);

                if (voiced[n] > kMarkerVoicedProbUpdate) {
                    auto trackedPeriod = 2.f / freq;
                    while (marker.period_ >
                           kMarkerSimilarPitchHighFactor * trackedPeriod) {
                        auto ratio  = marker.period_ / trackedPeriod;
                        auto iratio = std::round(ratio);
                        float newPeriod;
                        auto lastMarker = zffAnalysis_.getMarker();

                        // check if ratio is not integral
                        constexpr auto kRatioThresh = 0.2f;
                        if (std::abs(iratio - ratio) > kRatioThresh) {
                            newPeriod = lastMarker.period_;
                        } else {
                            newPeriod = marker.period_ / iratio;
                            newPeriod += kMarkerDetectedPitchRatio *
                                         (trackedPeriod - newPeriod);
                        }

                        DEBUGLOG << "oops \t" << trackedPeriod << " \t"
                                 << marker.period_ << " \t" << ratio << " \t"
                                 << iratio << " \t" << newPeriod << "\n";

                        assert(newPeriod > 0.f);
                        auto newPosition =
                            ctxt.modAdd(lastMarker.position_, newPeriod);
                        zffAnalysis_.addMarker({newPosition, newPeriod});
                        marker.period_ -= trackedPeriod;
                    }
                    if (marker.period_ >
                        kMarkerSimilarPitchLowFactor * trackedPeriod) {
                        // update period
                        marker.period_ += kMarkerDetectedPitchRatio *
                                          (trackedPeriod - marker.period_);
                        // update position
                        auto lastMarker = zffAnalysis_.getMarker();
                        marker.position_ =
                            ctxt.modAdd(lastMarker.position_, marker.period_);
                        assert(marker.period_ > 0.f);

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
            // smooth voiced enveloppe
            voicedSmooth_.process(dsp::Context(&voiced[n]), voicedSmoothState_);

            if (!voiced_ && voiced[n] > kPsolaHysteresisStart) {
                // set period of last marker
                auto period = 2.f / ctxtPitch.getInput();
                zffAnalysis_.getClosestMarker(ctxt).period_ = period;

                // start psola synthesis
                psola_.start(ctxt, zffAnalysis_);
                quantizer_.reset();
                voiced_ = true;
            } else if (voiced_ && voiced[n] < kPsolaHysteresisStop) {
                // stop psola synthesis
                psola_.stop();
                voiced_ = false;
            }
            if (voiced[n] > kGrainGenHysteresisStop) {
                grainGenerator_.stop();
            }
            if (voiced[n] < kGrainGenHysteresisStart) {
                grainGenerator_.start();
            }
            if (ctxtPitch.getInput() > 1e-6f) {
                auto inFreq = dsp::freq::NyquistRatio(ctxtPitch.getInput());
                auto pitch  = pitchRef_.toPitch(inFreq);
                auto quantizedPitch = quantizer_.process(pitch);
                auto pitchShift     = quantizedPitch - pitch;
                psola_.setFreq(inFreq.get());
                psola_.setShift(pitchShift.toRatio() * shift_);
                psola_.setPlayRate(formantShift_);
            }
            auto x = ctxt.getInput();

            // psola for voiced parts
            psola_.process(ctxt, synthDL_, zffAnalysis_);
            // grain generator for unvoiced parts
            grainGenerator_.process(dsp::BufferContext(&x, 1, ctxt.getBuffer()),
                                    synthDL_);
            // mix psola and grain generator
            x += (ctxt.getInput() - x) * voiced[n];
            ctxt.setOutput(x);

            ++n;
        };

        atomicShift_.store(psola_.getShift());

        count -= blockSize;
        buffer_.nextBlock(ctxt);
        mout += blockSize;
        pitch += blockSize;
        voiced += blockSize;
    }
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
