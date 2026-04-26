#include "PitchTracker.hpp"
#include <iostream>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T> static inline bool closeTo(T x, T y, T thresh)
{
    return std::abs(x - y) < thresh;
}

void PitchTracker::prepare(float sampleRate)
{
    auto freqScale = 2.f / sampleRate;
    freqScale_     = freqScale;

    decFactor_         = std::round(sampleRate / kTargetDownSampleRate);
    iDecFactor_        = static_cast<int>(decFactor_);
    minPeriodDetected_ = sampleRate / decFactor_ / kMaxDetectedFreq;

    // decimate nyquist filter
    // scipy.signal.cheby1(8,3,1,analog=True,output='sos')
    static constexpr decltype(nyquistFilter_)::SOS kSosAnalog = {
        {{{0.f, 0.f, 0.00783107f, 1.f, 0.21696145f, 0.05029392f},
          {0.f, 0.f, 1.f, 1.f, 0.18393103f, 0.32089197f},
          {0.f, 0.f, 1.f, 1.f, 0.12289879f, 0.7035754f},
          {0.f, 0.f, 1.f, 1.f, 0.04315631f, 0.97417345f}}}};

    nyquistFilter_.fromAnalog(kSosAnalog, freqScale * kMaxDetectedFreq);

    constexpr auto kAmpFollowerUp   = 0.0002f;
    constexpr auto kAmpFollowerDown = 0.01f;
    ampFollower_.setRate(0.1f, kAmpFollowerUp * sampleRate,
                         kAmpFollowerDown * sampleRate);

    constexpr auto kZcFollowerUp   = 0.002f;
    constexpr auto kZcFollowerDown = 0.005f;
    zeroCrossingFollower_.setRate(0.1f, kZcFollowerUp * sampleRate,
                                  kZcFollowerDown * sampleRate);

    setRange(range_);
}

void PitchTracker::setRange(Range range)
{
    switch (range) {
    default:
    case Range::kHigh:
        softLowpass_.setFreq(180.f * freqScale_);
        autoCorrelationThreshold_ = 0.84f;
        confirmPeriodMult_        = 2.1f;
        break;
    case Range::kLow:
        softLowpass_.setFreq(80.f * freqScale_);
        autoCorrelationThreshold_ = 0.97f;
        confirmPeriodMult_        = 4.f;
        break;
    }
}

dsp::interp::PolyPoint<float> PitchTracker::findNSDFPeriod()
{
    auto k = autoCorrelationThreshold_;

    auto nPeaks = nsdf_.updatePeaks(1e-8f);
    auto peaks  = nsdf_.getPeaks();

    // find max peak
    auto maxPeak = peaks[0];
    for (size_t i = 1; i < nPeaks; ++i) {
        if (peaks[i].y > maxPeak.y) {
            maxPeak = peaks[i];
        }
    }

    for (size_t i = 0; i < nPeaks; ++i) {

        // find best peak
        if (peaks[i].y >= k * maxPeak.y && peaks[i].x > minPeriodDetected_) {
            auto y      = peaks[i].y;
            auto period = peaks[i].x;

            // check harmonic jump
            if (voiced_ && downPeriod_ > 0.f) {
                auto ratio = downPeriod_ / period;

                HarmonicJump jump = kNoJump;
                float mult        = 1.f;
                auto newId        = i;
                if (closeTo(ratio, 2.f, 0.5f)) {
                    std::cout << "harmonic jump *2\n";
                    jump = kJumpOctaveUp;
                    mult = 2.f;
                    newId += 1;
                } else if (closeTo(ratio, 3.f, 0.5f)) {
                    std::cout << "harmonic jump *3\n";
                    jump = kJump3rdHarmonicUp;
                    mult = 3.f;
                    newId += 2;
                } else if (closeTo(ratio, 0.5f, 0.25f)) {
                    std::cout << "harmonic jump /2\n";
                    jump = kJumpOctaveDown;
                    mult = 0.5f;
                    newId -= 1;
                }

                if (jump != kNoJump) {
                    auto multPeriod = mult * period;
                    if (newId < nPeaks) {
                        auto newPeriod = peaks[newId].x;
                        if (closeTo((newPeriod - multPeriod) / period, 0.f,
                                    0.2f)) {
                            period = newPeriod;
                        } else {
                            period = multPeriod;
                        }
                    } else {
                        period = multPeriod;
                    }

                    ratio = downPeriod_ / period;
                    if (harmonicJumpType_ == jump) {
                        harmonicJumpCounter_ += 1;
                    } else {
                        harmonicJumpCounter_ = 1;
                    }
                    harmonicJumpType_ = jump;
                } else {
                    resetHarmonicJump();
                }

                auto thresh = static_cast<int>(period * 0.5);
                if (harmonicJumpCounter_ > thresh) {
                    period = peaks[i].x;
                    switch (harmonicJumpType_) {
                    case kJump3rdHarmonicUp:
                        std::cout << "change 3rd harmonic up\n";
                        downPeriod_ *= .33333333f;
                        break;
                    case kJumpOctaveDown:
                        std::cout << "change octave down\n";
                        downPeriod_ *= 2.f;
                        break;
                    case kJumpOctaveUp:
                        std::cout << "change octave up\n";
                        downPeriod_ *= .5f;
                        break;
                    default:
                        break;
                    }
                    resetHarmonicJump();
                }
            }
            return {period, y};
        }
    }

    constexpr auto kMinPeakY = -.5f;
    return {0, kMinPeakY};
}

float PitchTracker::smoothVoicedProbability(float p)
{
    p = std::clamp(p, 0.f, 1.f);
    p = hmm2_.process(p, 1.f - p);
    return p;
}

void PitchTracker::resetHarmonicJump()
{
    harmonicJumpCounter_ = 0;
    harmonicJumpType_    = kNoJump;
}

void PitchTracker::process(const float *const *__restrict in,
                           float *const *__restrict out, int count)
{
    const float *min = *in;
    float *pitchOut  = out[0];
    float *voicedOut = out[1];

    if (min != pitchOut) {
        for (int i = 0; i < count; ++i) pitchOut[i] = min[i];
    }
    if (min != voicedOut) {
        for (int i = 0; i < count; ++i) voicedOut[i] = min[i];
    }

    while (count) {
        auto blockSize = std::min(count, kMaxBlockSize);
        count -= blockSize;

        auto ctxtPitch  = dsp::BufferContext(pitchOut, blockSize, buffer_);
        auto ctxtVoiced = dsp::Context(voicedOut, blockSize);

        CTXTRUN(ctxtPitch)
        {
            nyquistFilter_.process(ctxtPitch, nyquistFilterState_);
        };

        // compute voiced probability
        CTXTRUN2(ctxtPitch, ctxtVoiced)
        {
            // follow input amplitude
            auto energy = ctxtVoiced.getInput();
            energy *= energy;
            ampFollower_.process(dsp::Context(&energy), ampState_);

            // follow low filtered amplitude
            auto lowEnergy = ctxtPitch.getInput();
            lowEnergy *= lowEnergy;
            ampFollower_.process(dsp::Context(&lowEnergy), ampFilteredState_);

            // follow zero crossing
            auto zc = zeroCrossing_.process(ctxtVoiced) == 0 ? 0.f : 1.f;
            zeroCrossingFollower_.process(dsp::Context(&zc),
                                          zeroCrossingFollowerState_);

            constexpr auto kEnergyProb0 = 0.0001f;
            constexpr auto kEnergyProb1 = 0.005f;
            auto hasEnergy =
                (energy - kEnergyProb0) / (kEnergyProb1 - kEnergyProb0);
            hasEnergy = std::clamp(hasEnergy, 0.f, 1.f);

            constexpr auto kEnergyRatioProb0 = 0.04f;
            constexpr auto kEnergyRatioProb1 = 0.20f;
            constexpr auto kEpsEnergy        = 0.000001f;
            auto energyRatio  = lowEnergy / (energy + kEpsEnergy);
            auto hasLowEnergy = (energyRatio - kEnergyRatioProb0) /
                                (kEnergyRatioProb1 - kEnergyRatioProb0);
            hasLowEnergy      = std::max(0.f, hasLowEnergy);

            constexpr auto kZCProb0 = 0.8f;
            constexpr auto kZCProb1 = 0.2f;
            auto hasZc              = (zc - kZCProb0) / (kZCProb1 - kZCProb0);
            hasZc                   = std::clamp(hasZc, 0.f, 1.f);

            auto voicedProb = hasEnergy * hasLowEnergy; //*hasZc;

            ctxtVoiced.setOutput(voicedProb);
        };

        // lowpass on signal to soften high harmonics
        CTXTRUN(ctxtPitch)
        {
            softLowpass_.process(ctxtPitch, softLowpassState_);
        };

        CTXTRUN2(ctxtPitch, ctxtVoiced)
        {
            if (voiced_) {
                nsdfTrack_.process(ctxtPitch, nsdfTrackDL_);

                if (nConfirmPitch_ > 0) {
                    // check if down period has shifted
                    auto periodDiff = period_ - nsdfTrack_.getPeriod();
                    if (std::abs(periodDiff) > 0.5f) {
                        if (periodDiff > 0.f)
                            nsdfTrack_.shiftUp(ctxtPitch, nsdfTrackDL_);
                        else
                            nsdfTrack_.shiftDown(ctxtPitch, nsdfTrackDL_);
                    }
                    // update voiced probability
                    auto voicedProb = ctxtVoiced.getInput();
                    voicedProb =
                        smoothVoicedProbability(voicedProb * downProb_);
                    ctxtVoiced.setOutput(voicedProb);
                    --nConfirmPitch_;

                    if (voicedProb < kVoicedProbabilityStopConfirm) {
                        voiced_ = false;
                        resetHarmonicJump();
                        std::cout << "stop confirm " << downProb_ << "\n";
                    }
                } else if (nConfirmPitch_ == 0) {
                    auto peak = nsdfTrack_.updateTrack(ctxtPitch, nsdfTrackDL_);
                    period_   = peak.x;

                    if (peak.y < kAutoCorrelationPeakStopTracking) {
                        std::cout << "stop track " << peak.y << "\n";
                        voiced_ = false;
                        // reset hmm
                        resetHarmonicJump();
                        downPeriod_ = period_ / decFactor_;
                        downProb_   = 0.6f;
                    } else {
                        auto voicedProb = ctxtVoiced.getInput();
                        voicedProb      = smoothVoicedProbability(voicedProb);
                        ctxtVoiced.setOutput(voicedProb);
                    }
                }
            }

            if (decimateId_ == 0) {
                // advance decimated buffer
                bufferDec_.nextId(1);
                auto ctxtDec = dsp::BufferContext(ctxtPitch, bufferDec_);
                nsdf_.process(ctxtDec, nsdfDL_);

                if (!voiced_ || nConfirmPitch_ > 0) {
                    auto peak = findNSDFPeriod();

                    // compute probability of peak
                    constexpr auto kPeakProb0 = -0.1f;
                    constexpr auto kPeakProb1 = 0.45f;
                    auto hasPeak =
                        (peak.y - kPeakProb0) / (kPeakProb1 - kPeakProb0);
                    hasPeak = std::max(0.f, hasPeak);

                    constexpr auto kDriftProb0 = 3.f;
                    constexpr auto kDriftProb1 = 0.2f;
                    auto drift = downPeriod_ > 0.f
                                     ? std::abs(peak.x - downPeriod_)
                                     : 0.f;
                    auto hasDrift =
                        (drift - kDriftProb0) / (kDriftProb1 - kDriftProb0);
                    hasDrift = std::clamp(hasDrift, 0.f, 1.f);

                    downPeriod_ = peak.x;
                    period_     = decFactor_ * downPeriod_;
                    downProb_   = hasPeak * hasDrift;
                }
            }

            if (!voiced_) {
                // update voiced probability
                auto voicedProb = ctxtVoiced.getInput();
                voicedProb = smoothVoicedProbability(voicedProb * downProb_);
                ctxtVoiced.setOutput(voicedProb);

                if (voicedProb > kVoicedProbabilityStartTracking) {
                    std::cout << "start " << period_ << "\n";
                    voiced_        = true;
                    nConfirmPitch_ = static_cast<int>(
                        std::ceil(confirmPeriodMult_ * period_));
                    nsdfTrack_.clear();
                    std::fill(std::begin(bufferData_), std::end(bufferData_),
                              0.f);
                    nsdfTrack_.setPeriod(period_);
                }
            }

            auto freq = period_ > 0.f ? 2.f / period_ : 0.f;
            ctxtPitch.setOutput(freq);

            decimateId_ = (decimateId_ + 1) % iDecFactor_;
        };

        buffer_.nextId(blockSize);

        pitchOut += blockSize;
        voicedOut += blockSize;
    }
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
