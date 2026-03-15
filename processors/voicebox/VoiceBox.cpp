#include "VoiceBox.h"

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

void VoiceBox::prepare(float sampleRate)
{
    pitchTracker_.prepare(sampleRate);

    auto freqScale = 2.f / sampleRate;
    dcblocker_.setFreq(50.f * freqScale);

    sampleRate_ = sampleRate;
    zffAnalysis_.setWindowSize(
        static_cast<int>(sampleRate * kAnalysisWindowSize));
}

void VoiceBox::update(float shift) { shift_ = shift; }

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

    auto pitch = pitchBuffer_;
    while (count) {
        auto blockSize = std::min(kMaxBlockSize, count);

        pitchTracker_.process(in, &pitch, blockSize);

        auto ctxt      = dsp::BufferContext(mout, blockSize, buffer_);
        auto ctxtPitch = dsp::Context(pitch, blockSize);

        CTXTRUN2(ctxt, ctxtPitch)
        {
            zffAnalysis_.process(ctxt, analysisDL_, [&, this](auto &marker) {
                auto freq = ctxtPitch.getInput();
                if (freq > 0.f) {
                    auto trackedPeriod = 2.f / freq;
                    while (marker.period_ > 1.05f * trackedPeriod) {
                        auto lastMarker = zffAnalysis_.getMarker();
                        auto newPosition =
                            ctxt.modAdd(lastMarker.position_, trackedPeriod);
                        zffAnalysis_.addMarker({newPosition, trackedPeriod});
                        marker.period_ -= trackedPeriod;
                    }
                    if (marker.period_ > 0.95f * trackedPeriod) {
                        auto lastMarker = zffAnalysis_.getMarker();
                        marker.period_ +=
                            0.9f * (trackedPeriod - marker.period_);
                        marker.position_ =
                            ctxt.modAdd(lastMarker.position_, marker.period_);
                        return true;
                    }
                    return false;
                }
                return true;
            });
        };

        CTXTRUN2(ctxt, ctxtPitch)
        {
            if (!voiced_ && ctxtPitch.getInput() > 0.f) {
                // if there is no marker find one
                if (zffAnalysis_.getNumMarkers() == 0) {
                    auto period = 2.f / ctxtPitch.getInput();
                    // find maximum peak
                    float peakVal{};
                    int peakId{};
                    for (int i = 1; i < static_cast<int>(period); ++i) {
                        auto val = analysisDL_.read(ctxt, i);
                        if (val > peakVal) {
                            peakId  = i;
                            peakVal = val;
                        }
                    }
                    zffAnalysis_.addMarker(
                        {ctxt.idDiff(static_cast<float>(peakId)), period});
                }

                // start psola synthesis
                psola_.start(ctxt, zffAnalysis_);
                voiced_ = true;
            } else if (voiced_ && ctxtPitch.getInput() == 0.f) {
                // stop psola synthesis
                psola_.stop(ctxt);
                zffAnalysis_.resetMarkers();
                voiced_ = false;
            }
            if (ctxtPitch.getInput() > 0.f) {
                // psola_.setFreq(ctxtPitch.getInput() * shift_);
                psola_.setFreq(atune(ctxtPitch.getInput() * shift_));
                // psola_.setFreq(402.f/24000.f*shift_);
            }
            psola_.process(ctxt, synthDL_, zffAnalysis_);
        };

        count -= blockSize;
        buffer_.nextBlock(ctxt);
        mout += blockSize;
        pitch += blockSize;
    }
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
