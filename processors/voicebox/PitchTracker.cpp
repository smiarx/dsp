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

    softLowpass_.setFreq(80.f * freqScale);

    // decimate nyquist filter
    // scipy.signal.cheby1(8,3,1,analog=True,output='sos')
    static constexpr decltype(nyquistFilter_)::SOS kSosAnalog = {
        {{{0.f, 0.f, 0.00783107f, 1.f, 0.21696145f, 0.05029392f},
          {0.f, 0.f, 1.f, 1.f, 0.18393103f, 0.32089197f},
          {0.f, 0.f, 1.f, 1.f, 0.12289879f, 0.7035754f},
          {0.f, 0.f, 1.f, 1.f, 0.04315631f, 0.97417345f}}}};

    nyquistFilter_.fromAnalog(kSosAnalog, freqScale * 1200);

    zeroCrossingFollower_.setRate(0.01f, 1.f / (0.002f * sampleRate));

    ampFollower_.setRate(0.01f, 1.f/(0.005f*sampleRate));
}

float PitchTracker::findNSDFPeriod(float zeroCrossing)
{
    constexpr auto kK = 0.98f;

    auto nPeaks = nsdf_.updatePeaks(1e-7f);
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
        if (peaks[i].y >= kK * maxPeak.y && peaks[i].x > 5.f) {
            auto y      = peaks[i].y;
            auto period = peaks[i].x;

            // set emission probability
            y = std::max(0.f, y - 0.1f);
            y /= (0.5f - 0.1f);

            bool hasJumped = false;
            if (voiced_ && downPeriod_ > 0.f) {
                auto ratio        = downPeriod_ / period;

                HarmonicJump jump = kNoJump;
                float mult        = 1.f;
                auto newId = i;
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

                if(jump != kNoJump)
                {
                    auto multPeriod = mult * period;
                    if(newId < nPeaks && newId >= 0)
                    {
                        auto newPeriod  = peaks[newId].x;
                        if (closeTo((newPeriod - multPeriod) / period, 0.f,
                                    0.2f)) {
                            period = newPeriod;
                        } else {
                            period = multPeriod;
                        }
                    }
                    else
                    {
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
                    hasJumped            = true;
                }
            }

            auto diffPeriod =
                std::max(0.f, 1.5f - std::abs(period - downPeriod_) * 0.1f);
            y *= diffPeriod;
            downPeriod_ = period;

            y *= std::max(0.f, 1.f - (zeroCrossing - 0.1f) / 0.3f);

            y = std::min(1.f, y);
            hmm2_.process(y, 1.f - y);
#ifndef NDEBUG
            //std::cout << "{" << period << ", " << peaks[i].y << "}\t"
            //          << zeroCrossing << " \t" << hmm2_.getState() << "\t";
            //std::cout << "[";
            //for (int j = 0; j < nPeaks; ++j) {
            //    std::cout << "{" << peaks[j].x << ", " << peaks[j].y << "}";
            //}
            //std::cout << "] \n";
#endif
            if (hmm2_.getState() > 0.5f) {
                if (hasJumped && voiced_) {
                    // init tracking
                    period_ = period * kDecFactor;
                    nConfirmPitch_ =
                        3 * static_cast<int>(std::ceil(period_));
                    nsdfTrack_.clear();
                    std::fill(std::begin(bufferData_), std::end(bufferData_),
                              0.f);
                    nsdfTrack_.setPeriod(period_);
                }
                return period;
            }
            return 0.f;
        }
    }

    return 0.f;
}

void PitchTracker::process(const float *const *__restrict in,
                           float *const *__restrict out, int count)
{
    const float *min = *in;
    float *mout      = *out;

    if (min != mout) {
        for (int i = 0; i < count; ++i) mout[i] = min[i];
    }

    while (count) {
        auto blockSize = std::min(count, kMaxBlockSize);
        count -= blockSize;

        auto ctxt   = dsp::BufferContext(mout, blockSize, buffer_);
        auto ctxtZc = dsp::Context(zeroCrossingBuffer_, blockSize);

        // compute zero crossing
        CTXTRUN2(ctxt, ctxtZc)
        {
            auto zc = zeroCrossing_.process(ctxt) == 0 ? 0.f : 1.f;
            zeroCrossingFollower_.process(dsp::Context(&zc),
                                          zeroCrossingFollowerState_);
            ctxtZc.setOutput(zc);
        };

        CTXTRUN(ctxt) { softLowpass_.process(ctxt, softLowpassState_); };

        CTXTRUN(ctxt) { nyquistFilter_.process(ctxt, nyquistFilterState_); };

        for (int n = 0; n < blockSize; ++n) {

            auto amp = ctxt.getInput();
            amp = std::abs(amp);
            ampFollower_.process(dsp::Context(&amp), ampFollowerState_);

            if (decimateId_ == 0) {
                bufferDec_.nextId(1);
                auto ctxtDec = dsp::BufferContext(ctxt, bufferDec_);
                nsdf_.process(ctxtDec, nsdfDL_);
                if (!voiced_) {
                    auto period = findNSDFPeriod(zeroCrossingBuffer_[n]);
                    period_ = kDecFactor * period;

                    if (period > 0) {
                        voiced_ = true;
                        nConfirmPitch_ =
                            3 * static_cast<int>(std::ceil(period_));
                        nsdfTrack_.clear();
                        std::fill(std::begin(bufferData_),
                                  std::end(bufferData_), 0.f);
                        nsdfTrack_.setPeriod(period_);
                    }
                } else if (nConfirmPitch_ > 0) {
                    auto period = findNSDFPeriod(zeroCrossingBuffer_[n]);

                    if (period > 0.f) {
                        downPeriod_ = period;
                        period_     = kDecFactor * period;
                    } else {
                        std::cout << "stop in test\n";
                        voiced_     = false;
                        decimateId_ = (decimateId_ + 1) % kDecFactor;
                    }
                }
            }

            if (voiced_) {
                nsdfTrack_.process(ctxt, nsdfTrackDL_);

                if (nConfirmPitch_ > 0) {
                    auto freq = 2.f / period_;
                    assert(!std::isnan(freq));
                    ctxt.setOutput(freq);
                    --nConfirmPitch_;
                } else {
                    auto peak = nsdfTrack_.updateTrack(ctxt, nsdfTrackDL_);
                    auto freq = 2.f / peak.x;
                    assert(!std::isnan(freq));
                    ctxt.setOutput(freq);

                    // probability
                    auto y  = peak.y;
                    y       = (y + 0.25f) * 2.f;
                    //auto zc = zeroCrossingBuffer_[n];
                    //y *= std::max(0.f, 1.f - (zc - 0.1f) / 0.4f);
                    y *= std::min(1.f, amp*1e2f);
                    if (y < 0.5f) {
#ifndef NDEBUG
                        std::cout << "stop " << peak.y << " amp:" << amp << "\n";
#endif
                        voiced_ = false;
                        // reset hmm
                        hmm2_.reset();
                        harmonicJumpCounter_ = 0;
                        harmonicJumpType_    = kNoJump;
                        downPeriod_          = 0.f;
                        ctxt.setOutput(0.f);
                    }
                }
            } else {
                ctxt.setOutput(0);
            }

            decimateId_ = (decimateId_ + 1) % kDecFactor;
            ctxt.next(1);
        }
        buffer_.nextId(blockSize);
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
