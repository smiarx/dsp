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
    dcblocker_.setFreq(50.f * freqScale);

    sampleRate_ = sampleRate;
    zffAnalysis_.setWindowSize(
        static_cast<int>(sampleRate * kAnalysisWindowSize));

    filter_.setWarpIn(barkScale(sampleRate));
    filter_.setWarpOut(barkScale(sampleRate));
}

void VoiceBox::update(float shift, float formantShift)
{
    if (shift != getPitchShift())
    {
        setPitchShift(shift);
    }
    if (formantShift != getFormantShift()) {
        setFormantShift(formantShift);
    }
}

void VoiceBox::setPitchShift(float shift)
{
    shift_ = shift;
}

void VoiceBox::setFormantShift(float formantShift)
{
        formantShift_ = formantShift;
        psola_.setPlayRate(formantShift);
        filter_.setWarpOut(barkScale(sampleRate_ / formantShift));
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

    auto pitch = pitchBuffer_;
    while (count) {
        auto blockSize = std::min(kMaxBlockSize, count);

        pitchTracker_.process(in, &pitch, blockSize);

        // set pitch atomic value to last pitch computed
        atomicFreq_.store(pitch[blockSize-1]*sampleRate_*0.5f);

        auto ctxt      = dsp::BufferContext(mout, blockSize, buffer_);
        auto ctxtPitch = dsp::Context(pitch, blockSize);

        CTXTRUN(ctxtPitch) {
            auto p = ctxtPitch.getInput();
        };

        CTXTRUN2(ctxt, ctxtPitch)
        {
            zffAnalysis_.process(ctxt, analysisDL_, [&, this](auto &marker) {
                // ctxt.setOutput(ctxt.getInput()-1.f);
                auto freq = ctxtPitch.getInput();
                assert(marker.period_ >= 0.f);

                if (freq > 0.f) {
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
                        assert(marker.period_ < 1000);

                        // ctxt.setOutput(ctxt.getInput()+2.f);
                        return true;
                    }
                    return false;
                }
                return true;
            });
        };

        CTXTRUN2(ctxt, ctxtPitch)
        {
            if (ctxtPitch.getInput() > 0.f) {
                psola_.setFreq(ctxtPitch.getInput() * shift_);
                // psola_.setFreq(atune(ctxtPitch.getInput() * shift_));
                // psola_.setFreq(400.f/24000.f*shift_);
            }
            if (!voiced_ && ctxtPitch.getInput() > 0.f) {
                //// if there is no marker find one
                //if (zffAnalysis_.getNumMarkers() == 0) {
                //    auto period = 2.f / ctxtPitch.getInput();
                //    // find maximum peak
                //    float peakVal{};
                //    int peakId{};
                //    for (int i = 1; i < static_cast<int>(period); ++i) {
                //        auto val = analysisDL_.read(ctxt, i);
                //        if (val > peakVal) {
                //            peakId  = i;
                //            peakVal = val;
                //        }
                //    }
                //    zffAnalysis_.addMarker(
                //        {ctxt.idDiff(static_cast<float>(peakId)), period});
                //}

                // set period of last marker
                auto period = 2.f / ctxtPitch.getInput();
                zffAnalysis_.getClosestMarker(ctxt).period_ = period;

                // start psola synthesis
                psola_.start(ctxt, zffAnalysis_);
                voiced_ = true;
            } else if (voiced_ && ctxtPitch.getInput() == 0.f) {
                // stop psola synthesis
                psola_.stop(ctxt);
                // zffAnalysis_.resetMarkers();
                voiced_ = false;
            }
            auto x = ctxt.getInput();
            synthDL_.write(ctxt, x);

            // x += noise_.process()*0.002f;
            // rls_.process(dsp::Context{&x}, warpAState1_, filter_);
            // auto cpyFilter = filter_;
            // cpyFilter.smoothCoeffs(0.92f);

            // cpyFilter.analyze(ctxt, warpAState_);
            // cpyFilter.compensateResidualAnalyze(ctxt, warpAState_);

            // cpyFilter.compensateResidualReconstruct(ctxt, warpRState_);
            // cpyFilter.reconstruct(ctxt, warpRState_);

            psola_.process(ctxt, synthDL_, zffAnalysis_);

            // ctxt.setOutput(noise_.process()*0.2f);
        };

        count -= blockSize;
        buffer_.nextBlock(ctxt);
        mout += blockSize;
        pitch += blockSize;
    }
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
