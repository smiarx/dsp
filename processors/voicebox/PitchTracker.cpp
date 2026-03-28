#include "PitchTracker.hpp"
#include <iostream>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

void PitchTracker::prepare(float sampleRate)
{
    auto freqScale = 2.f / sampleRate;

    constexpr auto kAmpDecay = 0.02f;
    ampFollower_.setRate(0.01f, 1.f / (kAmpDecay * sampleRate));

    softLowpass_.setFreq(80.f * freqScale);

    // decimate nyquist filter
    // scipy.signal.cheby1(8,3,1,analog=True,output='sos')
    static constexpr decltype(nyquistFilter_)::SOS kSosAnalog = {
        {{{0.f, 0.f, 0.00783107f, 1.f, 0.21696145f, 0.05029392f},
          {0.f, 0.f, 1.f, 1.f, 0.18393103f, 0.32089197f},
          {0.f, 0.f, 1.f, 1.f, 0.12289879f, 0.7035754f},
          {0.f, 0.f, 1.f, 1.f, 0.04315631f, 0.97417345f}}}};

    constexpr auto kNyquistFactor =
        0.9f; // we set the cutoff just above nyquist
    nyquistFilter_.fromAnalog(kSosAnalog,
                              kNyquistFactor / static_cast<float>(kDecFactor));

    zeroCrossingFollower_.setRate(0.01f, 1.f/(0.002f*sampleRate));
}

float PitchTracker::findNSDFPeriod(float zeroCrossing)
{
    constexpr auto kK             = 0.95f;
    constexpr auto kPeakThreshold = 0.45f;

    auto nPeaks = nsdf_.updatePeaks(1e-4f);
    auto peaks  = nsdf_.getPeaks();

    auto maxPeak = peaks[0];
    for (size_t i = 1; i < nPeaks; ++i) {
        if (peaks[i].y > maxPeak.y) {
            maxPeak = peaks[i];
        }
    }

    std::cout << "[ ";
    for(size_t j = 0; j < nPeaks; ++j)
    {
        std::cout << "{" << peaks[j].x << ", " << peaks[j].y << "}, ";
    }
    std::cout << "]\n";
    for (size_t i = 0; i < nPeaks; ++i) {
        if (peaks[i].y >= kK * maxPeak.y) {
            auto y = peaks[i].y;
            y=std::max(0.f, y-0.1f);
            y/=(0.5f-0.1f);

            static auto prevX = 0.f;
            auto x = peaks[i].x;
            auto diffx = std::max(0.f, 1.5f - std::abs(x-prevX)*0.25f);
            y *= diffx;
            prevX = x;

            y *= std::max(0.f,1.f-(zeroCrossing-0.03f)/0.3f);

            y = std::min(1.f,y);
            hmm2_.process(y,1.f-y);
            std::cout << "{" << peaks[i].x << ", " << peaks[i].y << "}\t" << zeroCrossing << "\t" << hmm2_.getState() << "\n";
            if(hmm2_.getState() > 0.5f) return peaks[i].x;
            //if (peaks[i].y > kPeakThreshold) return peaks[i].x;
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

    while(count)
    {
        auto blockSize = std::min(count, kMaxBlockSize);
        count -= blockSize;

        auto ctxt = dsp::Context(mout, blockSize);
        auto ctxtZc = dsp::Context(zeroCrossingBuffer_, blockSize);

        // compute zero crossing
        CTXTRUN2(ctxt, ctxtZc)
        {
            auto zc = zeroCrossing_.process(ctxt) == 0 ? 0.f : 1.f;
            zeroCrossingFollower_.process(dsp::Context(&zc), zeroCrossingFollowerState_);
            ctxtZc.setOutput(zc);
        };

        CTXTRUN(ctxt)
        {
            softLowpass_.process(ctxt, softLowpassState_);
            ctxt.setOutput(ctxt.getInput() * 10.f);
        };

        CTXTRUN(ctxt) { nyquistFilter_.process(ctxt, nyquistFilterState_); };

        while (blockSize > 0) {
            if (!voiced_) {
                auto nStart = (kDecFactor-decimateId_) % kDecFactor;
                auto ctxtDec = dsp::BufferContext(mout + nStart, blockSize, bufferDec_);
                int n;
                int decBlockSize = 0;
                for (n = nStart; n < ctxt.getBlockSize(); n += kDecFactor) {
                    nsdf_.process(ctxtDec, nsdfDL_);
                    auto period = findNSDFPeriod(zeroCrossingBuffer_[n]);

                    if (period > 0) {
                        auto freq = 2.f / static_cast<float>(kDecFactor) / period;
                        adaptiveNotchState_.setFreq(freq, 1e-6f);
                        voiced_            = true;
                        poleFollowerState_ = {};
                        //poleFollowerState_ = {1000};
                        ampFollowerState_ = {};
                        adaptiveNotchState_.setState(0.f);
                        nConfirmPitch_ = 2*static_cast<int>(std::ceil(period));
                        ++n;
                        break;
                    }

                    ctxtDec.next(kDecFactor, 1);
                    ++decBlockSize;
                }
                for (int i = 0; i < n; ++i) mout[i] = 0.f;
                decimateId_ = n > ctxt.getBlockSize() ? n % ctxt.getBlockSize() : n % kDecFactor;
                bufferDec_.nextId(decBlockSize);
                blockSize -= n;
                mout += n;
            }

            else if (voiced_) {
                auto ctxt = dsp::Context(mout, blockSize);
                int n;
                
                if(nConfirmPitch_ > 0)
                {
                    for (n = 0; n < ctxt.getBlockSize() && nConfirmPitch_ > 0; n += 1, ctxt.next()) {
                        if(decimateId_ == 0)
                        {
                            // set correct buffer id
                            bufferDec_.nextId(1);
                            nsdf_.process(dsp::BufferContext(ctxt, bufferDec_), nsdfDL_);
                            auto period = findNSDFPeriod(zeroCrossingBuffer_[n]);
                            
                            if(period > 0.f)
                            {
                                auto freq = 2.f / static_cast<float>(kDecFactor) / period;
                                adaptiveNotchState_.setFreq(freq);
                            }
                            else
                            {
                                voiced_ = false;
                                break;
                            }
                        }

                        auto amp = ctxt.getInput();
                        amp *= amp;

                        auto x = ctxt.getInput();
                        auto a = adaptiveNotchState_.getCoeff();
                        notchState_.setCoeff(a*a-2.f);
                        notch_.process(ctxt, notchState_);

                        adaptiveNotch_.process(ctxt, adaptiveNotchState_);
                        auto freq = adaptiveNotchState_.getFreq();

                        //auto a = adaptiveNotchState_.getCoeff();
                        //notchState_.setCoeff(a);
                        //auto pole = x;
                        //notch_.process(dsp::Context(&pole), notchState_);

                        auto pole = ctxt.getInput();


                        //auto pole = std::abs(adaptiveNotchState_.getState()[0]);
                        //pole *= 0.002f;

                        pole *= pole;

                        ampFollower_.process(dsp::Context(&amp), ampFollowerState_);
                        ampFollower_.process(dsp::Context(&pole), poleFollowerState_);

                        ctxt.setOutput(freq);
                        --nConfirmPitch_;
                        decimateId_ = (decimateId_ + 1) % kDecFactor;
                    }
                }
                else
                {
                    for (n = 0; n < ctxt.getBlockSize(); n += 1, ctxt.next()) {
                        auto amp = ctxt.getInput();
                        amp      = std::abs(amp);

                        auto x = ctxt.getInput();
                        auto a = adaptiveNotchState_.getCoeff();
                        notchState_.setCoeff(a*a-2.f);
                        notch_.process(ctxt, notchState_);
                        adaptiveNotch_.process(ctxt, adaptiveNotchState_);
                        auto freq = adaptiveNotchState_.getFreq();
                        // std::cout << 24000*freq << "\n";

                        //auto a = adaptiveNotchState_.getCoeff();
                        //notchState_.setCoeff(a);
                        //auto pole = x;
                        //notch_.process(dsp::Context(&pole), notchState_);

                        auto pole = ctxt.getInput();

                        //auto pole = std::abs(adaptiveNotchState_.getState()[0]);
                        //pole *= 0.002f;

                        pole = std::abs(pole);

                        ampFollower_.process(dsp::Context(&amp), ampFollowerState_);
                        ampFollower_.process(dsp::Context(&pole), poleFollowerState_);

                        if (pole > 0.99f * amp || amp < 1e-4f) {
                        //if (pole < amp || amp < 0.00001f) {
                            voiced_ = false;
                            ctxt.setOutput(0.f);
                            nsdf_.clear();
                            std::fill(std::begin(bufferDecData_),
                                      std::end(bufferDecData_), 0.f);
                            n++;
                            break;
                        }

                        ctxt.setOutput(freq);
                    };
                }
                decimateId_ = n % kDecFactor;
                blockSize -= n;
                mout += n;
            }
        }
    }
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
