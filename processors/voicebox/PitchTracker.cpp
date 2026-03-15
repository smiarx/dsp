#include "PitchTracker.hpp"

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
{

void PitchTracker::prepare(float sampleRate)
{
    auto freqScale = 2.f / sampleRate;

    constexpr auto kAmpDecay = 0.005f;
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
}

float PitchTracker::findNSDFPeriod()
{
    constexpr auto kK             = 0.92f;
    constexpr auto kPeakThreshold = 0.48f;

    auto nPeaks = nsdf_.updatePeaks();
    auto peaks  = nsdf_.getPeaks();

    auto maxPeak = peaks[0];
    for (size_t i = 1; i < nPeaks; ++i) {
        if (peaks[i].y > maxPeak.y) {
            maxPeak = peaks[i];
        }
    }

    for (size_t i = 0; i < nPeaks; ++i) {
        if (peaks[i].y >= kK * maxPeak.y) {
            if (peaks[i].y > kPeakThreshold) return peaks[i].x;
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

    auto ctxt = dsp::Context(mout, count);

    CTXTRUN(ctxt)
    {
        softLowpass_.process(ctxt, softLowpassState_);
        ctxt.setOutput(ctxt.getInput() * 10.f);
    };

    CTXTRUN(ctxt) { nyquistFilter_.process(ctxt, nyquistFilterState_); };

    while (count > 0) {
        if (!voiced_) {
            auto ctxtDec = dsp::BufferContext(mout, count, bufferDec_);
            int n;
            int decBlockSize = 0;
            for (n = decimateId_; n < ctxt.getBlockSize(); n += kDecFactor) {
                nsdf_.process(ctxtDec, nsdfDL_);
                auto period = findNSDFPeriod();
                float freq  = 0.f;
                if (period > 0) {
                    freq = 2.f / static_cast<float>(kDecFactor) / period;
                    adaptiveNotchState_.setFreq(freq);
                    voiced_            = true;
                    poleFollowerState_ = 10000;
                    break;
                }

                ctxtDec.next(kDecFactor, 1);
                ++decBlockSize;
            }
            for (int i = 0; i < n; ++i) mout[i] = 0.f;
            decimateId_ = n % ctxt.getBlockSize();
            bufferDec_.nextId(decBlockSize);
            count -= n;
            mout += n;
        }

        if (voiced_) {
            auto ctxt = dsp::Context(mout, count);
            int n;
            for (n = 0; n < ctxt.getBlockSize(); n += 1, ctxt.next()) {
                auto amp = ctxt.getInput();
                amp      = amp * amp;

                adaptiveNotch_.process(ctxt, adaptiveNotchState_);
                auto freq = adaptiveNotchState_.getFreq();
                // std::cout << 24000*freq << "\n";

                auto pole = std::abs(adaptiveNotchState_.getPole());
                pole *= 0.002f;
                pole *= pole;

                ampFollower_.process(dsp::Context(&amp), ampFollowerState_);
                ampFollower_.process(dsp::Context(&pole), poleFollowerState_);

                if (pole < amp || amp < 0.0001f) {
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
            decimateId_ = n % ctxt.getBlockSize();
            count -= n;
            mout += n;
        }
    }
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
