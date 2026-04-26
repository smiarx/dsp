#include "VoiceBox.h"
#include <iostream>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

void VoiceBox::prepare(float sampleRate)
{
    pitchTracker_.prepare(sampleRate);

    auto freqScale = 2.f / sampleRate;

    sampleRate_ = sampleRate;
    freqScale_  = freqScale;

    zffAnalysis_.setWindowSize(
        static_cast<int>(sampleRate * kAnalysisWindowSize));

    /* signal.cheby2(8,80,1,"high", analog=True, output="sos") */
    static constexpr decltype(dcblocker_)::SOS kSosAnalog = {
        {{{1.f, 0.f, 0.03806023f, 1.f, 2.15704486f, 1.24729461f},
          {1.f, 0.f, 0.30865828f, 1.f, 1.82865434f, 1.51789266f},
          {1.f, 0.f, 0.69134172f, 1.f, 1.82865434f, 1.51789266f},
          {1.f, 0.f, 0.96193977f, 1.f, 1.82865434f, 1.51789266f}}}};
    dcblocker_.fromAnalog(kSosAnalog, 40.f * freqScale);

    grainGenerator_.setFreq(100.f * freqScale);
    grainGenerator_.start();

    voicedSmooth_.setRate(0.01f, 1.f / (sampleRate * 0.05f));

    setCorrectSpeed(correctSpeed_);
    setA4Frequency(a4Frequency_);
    setScale(scale_);
}

void VoiceBox::update(float shift, float formantShift, float correctSpeed,
                      float a4Frequency)
{
    if (shift != getPitchShift()) {
        setPitchShift(shift);
    }
    if (formantShift != getFormantShift()) {
        setFormantShift(formantShift);
    }
    if (correctSpeed != getCorrectSpeed()) {
        setCorrectSpeed(correctSpeed);
    }
    if (a4Frequency != getA4Frequency()) {
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

void VoiceBox::setCorrectSpeed(float time)
{
    std::cout << time << "\n";
    correctSpeed_ = time;
    correctionFilter_.setRate(0.01f, 1.f / (sampleRate_ * time));
}

void VoiceBox::setA4Frequency(float freq)
{
    a4Frequency_ = freq;
    pitchBias_   = std::log2(freq * freqScale_);
}

void VoiceBox::setScale(Scale scale)
{
    scale_         = scale;
    nScalePitches_ = 0;
    for (int n = 0; n < kNumNotes; ++n) {
        auto note = static_cast<Note>(n);
        if (scale_ & note) {
            scalePitches_[nScalePitches_] =
                static_cast<float>(n) / static_cast<float>(kNumNotes);
            ++nScalePitches_;
        }
    }
    scalePitches_[nScalePitches_] = scalePitches_[0] + 1.f;
    ++nScalePitches_;

    for (int i = 0; i < nScalePitches_; ++i) {
        std::cout << scalePitches_[i] << " - ";
    }
    std::cout << "\n";
}

void VoiceBox::setRange(PitchTracker::Range range)
{
    pitchTracker_.setRange(range);
}

float VoiceBox::computePitchShift(float inFreq)
{
    auto pitch         = std::log2(inFreq) - pitchBias_;
    auto octave        = std::floor(pitch);
    auto pitchInOctave = pitch - octave;

    float minDiff     = std::abs(scalePitches_[0] - pitchInOctave);
    float chosenPitch = scalePitches_[0];
    for (size_t i = 1; i < nScalePitches_; ++i) {
        auto diff = std::abs(scalePitches_[i] - pitchInOctave);
        if (diff < minDiff) {
            minDiff     = diff;
            chosenPitch = scalePitches_[i];
        }
    }

    chosenPitch += octave;

    // smooth filter pitch
    correctionFilter_.process(dsp::Context(&chosenPitch),
                              correctionFilterState_);

    auto shift = std::exp2(chosenPitch - pitch);
    return shift;
};

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
                        std::cout << "oops \t" << trackedPeriod << " \t"
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

            if (ctxtPitch.getInput() > 1e-6f) {
                auto inFreq = ctxtPitch.getInput();
                auto shift  = computePitchShift(inFreq);
                psola_.setFreq(inFreq);
                psola_.setShift(shift * shift_);
                psola_.setPlayRate(shift * formantShift_);
            }
            if (!voiced_ && voiced[n] > kPsolaHysteresisStart) {
                // set period of last marker
                auto period = 2.f / ctxtPitch.getInput();
                zffAnalysis_.getClosestMarker(ctxt).period_ = period;

                // start psola synthesis
                psola_.start(ctxt, zffAnalysis_);
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

        count -= blockSize;
        buffer_.nextBlock(ctxt);
        mout += blockSize;
        pitch += blockSize;
        voiced += blockSize;
    }
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
