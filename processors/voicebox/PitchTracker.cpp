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

    softLowpass_.setFreq(180.f * freqScale);

    // decimate nyquist filter
    // scipy.signal.cheby1(8,3,1,analog=True,output='sos')
    static constexpr decltype(nyquistFilter_)::SOS kSosAnalog = {
        {{{0.f, 0.f, 0.00783107f, 1.f, 0.21696145f, 0.05029392f},
          {0.f, 0.f, 1.f, 1.f, 0.18393103f, 0.32089197f},
          {0.f, 0.f, 1.f, 1.f, 0.12289879f, 0.7035754f},
          {0.f, 0.f, 1.f, 1.f, 0.04315631f, 0.97417345f}}}};

    nyquistFilter_.fromAnalog(kSosAnalog, freqScale * 1200);

    constexpr auto kAmpFollowerUp   = 0.0002f;
    constexpr auto kAmpFollowerDown = 0.01f;
    ampFollower_.setRate(0.1f, kAmpFollowerUp * sampleRate,
                         kAmpFollowerDown * sampleRate);

    constexpr auto kZcFollowerUp   = 0.002f;
    constexpr auto kZcFollowerDown = 0.005f;
    zeroCrossingFollower_.setRate(0.1f, kZcFollowerUp * sampleRate,
                                  kZcFollowerDown * sampleRate);

    smoothFreq_.setFreq(60 * freqScale);
}

void PitchTracker::setRange(Range range)
{
    switch (range) {
    default:
    case Range::kHigh:
        softLowpass_.setFreq(180.f * freqScale_);
        autoCorrecationThreshold_ = 0.84f;
        break;
    case Range::kLow:
        softLowpass_.setFreq(80.f * freqScale_);
        autoCorrecationThreshold_ = 0.97f;
        break;
    }
}

dsp::interp::PolyPoint<float> PitchTracker::findNSDFPeriod()
{
    auto k = autoCorrecationThreshold_;

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
        // TODO change fixed value
        if (peaks[i].y >= k * maxPeak.y && peaks[i].x > 5.f) {
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
                    if (newId < nPeaks && newId >= 0) {
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
                    harmonicJumpCounter_ = 0;
                    harmonicJumpType_    = kNoJump;
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
                    harmonicJumpCounter_ = 0;
                    harmonicJumpType_    = kNoJump;
                }
            }
            return {period, y};
        }
    }

    return {};
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
                    // update voiced probability
                    auto voicedProb = ctxtVoiced.getInput();
                    voicedProb *= downProb_;
                    voicedProb = std::clamp(voicedProb, 0.f, 1.f);
                    voicedProb *= 0.99f;
                    voicedProb = hmm2_.process(voicedProb, 1.f - voicedProb);
                    ctxtVoiced.setOutput(voicedProb);
                    --nConfirmPitch_;

                    if (voicedProb < 0.2f) {
                        voiced_              = false;
                        harmonicJumpCounter_ = 0;
                        harmonicJumpType_    = kNoJump;
                        std::cout << "stop confirm " << downProb_ << "\n";
                    }
                } else if (nConfirmPitch_ == 0) {
                    auto peak = nsdfTrack_.updateTrack(ctxtPitch, nsdfTrackDL_);
                    period_   = peak.x;

                    if (peak.y < 0.1f) {
                        std::cout << "stop track " << peak.y << "\n";
                        voiced_ = false;
                        // reset hmm
                        harmonicJumpCounter_ = 0;
                        harmonicJumpType_    = kNoJump;
                        downPeriod_          = period_ / kDecFactor;
                        downProb_            = 0.6f;
                    } else {
                        auto voicedProb = ctxtVoiced.getInput();
                        voicedProb      = std::clamp(voicedProb, 0.f, 1.f);
                        voicedProb *= 0.99f;
                        voicedProb =
                            hmm2_.process(voicedProb, 1.f - voicedProb);
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

                    if (peak.x < 4.f) {
                        downProb_ = 0.f;
                    } else {
                        downPeriod_ = peak.x;
                        period_     = kDecFactor * downPeriod_;
                        downProb_   = hasPeak * hasDrift;
                    }
                }
            }

            if (!voiced_) {
                // update voiced probability
                auto voicedProb = ctxtVoiced.getInput();
                voicedProb *= downProb_;
                voicedProb = std::clamp(voicedProb, 0.f, 1.f);
                voicedProb *= 0.99f;
                voicedProb = hmm2_.process(voicedProb, 1.f - voicedProb);
                ctxtVoiced.setOutput(voicedProb);

                if (voicedProb > 0.8f) {
                    std::cout << "start " << period_ << "\n";
                    voiced_        = true;
                    nConfirmPitch_ = 3 * static_cast<int>(std::ceil(period_));
                    nsdfTrack_.clear();
                    std::fill(std::begin(bufferData_), std::end(bufferData_),
                              0.f);
                    nsdfTrack_.setPeriod(period_);
                }
            }

            else {
            }

            auto freq = period_ > 0.f ? 2.f / period_ : 0.f;
            // TODO smooth req but reset when entering a new voiced sequence
            // smoothFreq_.process(dsp::Context(&freq), smoothFreqState_);
            ctxtPitch.setOutput(freq);

            decimateId_ = (decimateId_ + 1) % kDecFactor;
        };

        buffer_.nextId(blockSize);

        pitchOut += blockSize;
        voicedOut += blockSize;
    }
}

// void PitchTracker::process(const float *const *__restrict in,
//                            float *const *__restrict out, int count)
//{
//     const float *min = *in;
//     float *mout      = *out;
//
//     if (min != mout) {
//         for (int i = 0; i < count; ++i) mout[i] = min[i];
//     }
//
//     while(count)
//     {
//         auto blockSize = std::min(count, kMaxBlockSize);
//         count -= blockSize;
//
//         auto ctxt = dsp::Context(mout, blockSize);
//         auto ctxtZc = dsp::Context(zeroCrossingBuffer_, blockSize);
//
//         // compute zero crossing
//         CTXTRUN2(ctxt, ctxtZc)
//         {
//             auto zc = zeroCrossing_.process(ctxt) == 0 ? 0.f : 1.f;
//             zeroCrossingFollower_.process(dsp::Context(&zc),
//             zeroCrossingFollowerState_); ctxtZc.setOutput(zc);
//         };
//
//         CTXTRUN(ctxt)
//         {
//             softLowpass_.process(ctxt, softLowpassState_);
//             //ctxt.setOutput(ctxt.getInput() * 10.f);
//         };
//
//         CTXTRUN(ctxt) { nyquistFilter_.process(ctxt, nyquistFilterState_); };
//
//         while (blockSize > 0) {
//             if (!voiced_) {
//                 auto nStart = (kDecFactor-decimateId_) % kDecFactor;
//                 auto ctxtDec = dsp::BufferContext(mout + nStart, blockSize,
//                 bufferDec_); int n; int decBlockSize = 0; for (n = nStart; n
//                 < ctxt.getBlockSize(); n += kDecFactor) {
//                     nsdf_.process(ctxtDec, nsdfDL_);
//                     auto period = findNSDFPeriod(zeroCrossingBuffer_[n]);
//
//                     if (period > 0) {
//                         auto freq = 2.f / static_cast<float>(kDecFactor) /
//                         period; adaptiveNotchState_.setFreq(freq, 1e-4f);
//                         voiced_            = true;
//                         poleFollowerState_ = {};
//                         //poleFollowerState_ = {1000};
//                         ampFollowerState_ = {};
//                         adaptiveNotchState_.setState(0.f);
//                         nConfirmPitch_ =
//                         2*kDecFactor*static_cast<int>(std::ceil(period));
//                         ++n;
//                         break;
//                     }
//
//                     ctxtDec.next(kDecFactor, 1);
//                     ++decBlockSize;
//                 }
//                 n = std::min(n, blockSize);
//                 for (int i = 0; i < n; ++i) mout[i] = 0.f;
//                 decimateId_ = n > ctxt.getBlockSize() ? n %
//                 ctxt.getBlockSize() : n % kDecFactor;
//                 bufferDec_.nextId(decBlockSize);
//                 blockSize -= n;
//                 mout += n;
//             }
//
//             else if (voiced_) {
//                 auto ctxt = dsp::Context(mout, blockSize);
//                 int n;
//
//                 if(nConfirmPitch_ > 0)
//                 {
//                     for (n = 0; n < ctxt.getBlockSize() && nConfirmPitch_ >
//                     0; n += 1, ctxt.next()) {
//                         if(decimateId_ == 0)
//                         {
//                             // set correct buffer id
//                             bufferDec_.nextId(1);
//                             nsdf_.process(dsp::BufferContext(ctxt,
//                             bufferDec_), nsdfDL_); auto period =
//                             findNSDFPeriod(zeroCrossingBuffer_[n]);
//
//                             if(period > 0.f)
//                             {
//                                 auto freq = 2.f /
//                                 static_cast<float>(kDecFactor) / period;
//                                 adaptiveNotchState_.setFreq(freq);
//                             }
//                             else
//                             {
//                                 voiced_ = false;
//                                 break;
//                             }
//                         }
//
//                         auto a = adaptiveNotchState_.getCoeff();
//
//                         notchState_.setCoeff(a*a-2.f);
//                         notch_.process(ctxt, notchState_);
//                         notchState2_.setCoeff(a*(a*a-3.f));
//                         notch_.process(ctxt, notchState2_);
//
//                         constexpr auto r = 0.97f;
//                         {
//                         auto s = ctxt.getInput() + a*r*lp_[0] - (r*r)*lp_[1];
//                         auto y = (1.f-r)*(s + lp_[1]);
//                         lp_[1] = lp_[0];
//                         lp_[0] = s;
//                         ctxt.setOutput(y);
//                         }
//                         {
//                         auto s = ctxt.getInput() + a*r*lp2_[0] -
//                         (r*r)*lp2_[1]; auto y = (1.f-r)*(s + lp2_[1]);
//                         lp2_[1] = lp2_[0];
//                         lp2_[0] = s;
//                         ctxt.setOutput(y);
//                         }
//
//                         auto amp = ctxt.getInput();
//                         amp *= amp;
//
//                         auto x = ctxt.getInput();
//
//                         auto prevState = adaptiveNotchState_;
//                         adaptiveNotch_.process(ctxt, adaptiveNotchState_);
//                         adaptiveNotchState_.setCoeff(prevState.getCoeff());
//                         auto freq = adaptiveNotchState_.getFreq();
//
//                         //auto a = adaptiveNotchState_.getCoeff();
//                         //notchState_.setCoeff(a);
//                         //auto pole = x;
//                         //notch_.process(dsp::Context(&pole), notchState_);
//                         //
//                         //
//
//                         auto pole = ctxt.getInput();
//
//
//                         //auto pole =
//                         std::abs(adaptiveNotchState_.getState()[0]);
//                         //pole *= 0.002f;
//
//                         pole *= pole;
//
//                         ampFollower_.process(dsp::Context(&amp),
//                         ampFollowerState_);
//                         ampFollower_.process(dsp::Context(&pole),
//                         poleFollowerState_);
//
//                         ctxt.setOutput(freq);
//                         --nConfirmPitch_;
//                         decimateId_ = (decimateId_ + 1) % kDecFactor;
//                     }
//                 }
//                 else
//                 {
//                     for (n = 0; n < ctxt.getBlockSize(); n += 1, ctxt.next())
//                     {
//
//                         auto a = adaptiveNotchState_.getCoeff();
//                         notchState_.setCoeff(a*a-2.f);
//                         notch_.process(ctxt, notchState_);
//                         notchState2_.setCoeff(a*(a*a-3.f));
//                         notch_.process(ctxt, notchState2_);
//
//                         constexpr auto r = 0.97f;
//                         {
//                         auto s = ctxt.getInput() + a*r*lp_[0] - (r*r)*lp_[1];
//                         auto y = (1.f-r)*(s + lp_[1]);
//                         lp_[1] = lp_[0];
//                         lp_[0] = s;
//                         ctxt.setOutput(y);
//                         }
//                         {
//                         auto s = ctxt.getInput() + a*r*lp2_[0] -
//                         (r*r)*lp2_[1]; auto y = (1.f-r)*(s + lp2_[1]);
//                         lp2_[1] = lp2_[0];
//                         lp2_[0] = s;
//                         ctxt.setOutput(y);
//                         }
//
//                         auto amp = ctxt.getInput();
//                         amp      = std::abs(amp);
//
//                         auto x = ctxt.getInput();
//                         adaptiveNotch_.process(ctxt, adaptiveNotchState_);
//                         auto freq = adaptiveNotchState_.getFreq();
//                         // std::cout << 24000*freq << "\n";
//
//                         //auto a = adaptiveNotchState_.getCoeff();
//                         //notchState_.setCoeff(a);
//                         //auto pole = x;
//                         //notch_.process(dsp::Context(&pole), notchState_);
//
//                         auto pole = ctxt.getInput();
//
//
//                         //auto pole =
//                         std::abs(adaptiveNotchState_.getState()[0]);
//                         //pole *= 0.002f;
//
//                         pole = std::abs(pole);
//
//                         ampFollower_.process(dsp::Context(&amp),
//                         ampFollowerState_);
//                         ampFollower_.process(dsp::Context(&pole),
//                         poleFollowerState_);
//
//                         {
//                         auto zc = zeroCrossingBuffer_[n];
//                         auto y =
//                         std::max(0.f,1.f-(ctxtZc.getInput()-0.03f)/0.44f);
//                         hmm2_.process(y,1.f-y);
//                         }
//
//                         //if (hmm2_.getState() < 0.5f || amp < 1e-4f) {
//                         //if (pole < amp || amp < 0.00001f) {
//                         //if(pole > 0.9f*amp || amp < 0.0001f) {
//                         if(pole > 0.5f*amp || amp < 0.001f) {
//                             voiced_ = false;
//                             ctxt.setOutput(0.f);
//                             nsdf_.clear();
//                             std::fill(std::begin(bufferDecData_),
//                                       std::end(bufferDecData_), 0.f);
//                             n++;
//                             break;
//                         }
//
//                         ctxt.setOutput(freq);
//                         //ctxt.setOutput(x);
//                     };
//                 }
//                 decimateId_ = n % kDecFactor;
//                 blockSize -= n;
//                 mout += n;
//             }
//         }
//     }
// }

// void PitchTracker::prepare(float sampleRate)
//{
//     auto freqScale = 2.f / sampleRate;
//     softLowpass_.setFreq(400.f * freqScale);
//
//     dcRemove_.butterworthHP(60*freqScale);
//
//     // decimate nyquist filter
//     // scipy.signal.cheby1(8,3,1,analog=True,output='sos')
//     static constexpr decltype(nyquistFilter_)::SOS kSosAnalog = {
//         {{{0.f, 0.f, 0.00783107f, 1.f, 0.21696145f, 0.05029392f},
//           {0.f, 0.f, 1.f, 1.f, 0.18393103f, 0.32089197f},
//           {0.f, 0.f, 1.f, 1.f, 0.12289879f, 0.7035754f},
//           {0.f, 0.f, 1.f, 1.f, 0.04315631f, 0.97417345f}}}};
//
//     nyquistFilter_.fromAnalog(kSosAnalog,
//                               1000*freqScale);
//
//     zeroCrossingFollower_.setRate(0.01f, 1.f/(0.01f*sampleRate));
// }
//
// void PitchTracker::process(const float *const *__restrict in,
//                            float *const *__restrict out, int count)
//{
//     const float *min = *in;
//     float *mout      = *out;
//
//     if (min != mout) {
//         for (int i = 0; i < count; ++i) mout[i] = min[i];
//     }
//
//     while(count)
//     {
//         auto blockSize = std::min(count, kMaxBlockSize);
//         count -= blockSize;
//
//         auto ctxt = dsp::Context(mout, blockSize);
//         auto ctxtZc = dsp::Context(zeroCrossingBuffer_, blockSize);
//
//         // compute zero crossing
//         CTXTRUN2(ctxt, ctxtZc)
//         {
//             auto zc = zeroCrossing_.process(ctxt) == 0 ? 0.f : 1.f;
//             zeroCrossingFollower_.process(dsp::Context(&zc),
//             zeroCrossingFollowerState_); ctxtZc.setOutput(zc);
//         };
//
//         CTXTRUN(ctxt)
//         {
//             softLowpass_.process(ctxt, softLowpassState_);
//             //ctxt.setOutput(ctxt.getInput() * 10.f);
//
//             dcRemove_.process(ctxt, dcRemoveState_);
//
//             nyquistFilter_.process(ctxt, nyquistFilterState_);
//         };
//
//
//         CTXTRUN2(ctxt,ctxtZc)
//         {
//             auto a = adaptiveNotchState_.getCoeff();
//
//             auto y = std::max(0.f,1.f-(ctxtZc.getInput()-0.03f)/0.44f);
//             y *= std::max(0.f, 1.f-std::abs(((a-prevA_)*1000.f)));
//             hmm2_.process(y,1.f-y);
//             prevA_ = a;
//
//             {
//                 notchState_.setCoeff(a*a-2.f);
//                 notch_.process(ctxt, notchState_);
//
//                 constexpr auto r = 0.99f;
//                 auto s = ctxt.getInput() + a*r*lp_[0] - (r*r)*lp_[1];
//                 auto y = (1.f-r)*(1.f-r)*(s + lp_[0]+lp_[0] - lp_[1]);
//                 lp_[1] = lp_[0];
//                 lp_[0] = s;
//                 ctxt.setOutput(y);
//
//                 adaptiveNotch_.process(ctxt, adaptiveNotchState_);
//
//                 adaptiveNotch_.process(ctxt, adaptiveNotchState2_);
//             }
//                 auto freq = adaptiveNotchState_.getFreq();
//                 //freq = std::min(freq, adaptiveNotchState2_.getFreq());
//
//
//             if(hmm2_.getState() > 0.5f)
//             {
//                 ctxt.setOutput(freq);
//
//             }
//             else
//             {
//                 ctxt.setOutput(0.f);
//             }
//             //ctxt.setOutput(hmm2_.getState());
//             //ctxt.setOutput(ctxtZc.getInput());
//         };
//
//         mout += blockSize;
//     }
// }

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
