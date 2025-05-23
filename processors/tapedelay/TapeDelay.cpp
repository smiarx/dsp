#include "TapeDelay.h"
#include "dsp/Buffer.h"
#include "dsp/Context.h"
#include "dsp/FastMath.h"
#include "dsp/Signal.h"
#include "dsp/Utils.h"
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace processors
{
TapeDelay::FadeLut TapeDelay::fadeLut;

void TapeDelay::update(float delay, float feedback, float cutlowpass,
                       float cuthighpass, float saturation, float drift,
                       Mode mode, float drywet, int blockSize)
{
    if (delay != delay_) {
        setDelay(delay, blockSize);
    }
    if (drift != drift_) {
        setDrift(drift, blockSize);
    }
    if (cutlowpass != cutlowpass_) {
        setCutLowPass(cutlowpass, blockSize);
    }
    if (cuthighpass != cuthighpass_) {
        setCutHiPass(cuthighpass, blockSize);
    }
    if (mode != mode_) {
        setMode(mode, blockSize);
    }
    setSaturation(saturation, blockSize);
    setFeedback(feedback, blockSize);
    setDryWet(drywet, blockSize);
}

void TapeDelay::setDelay(float delay, int blockSize)
{
    delay_ = delay;

    delay *= sampleRate_;
    // limits;
    delay = std::max(delay, 0.f);
    if (mode_ == kReverse) {
        delay = std::min(delay, kMaxDelay * kDefaultSampleRate /
                                    kReverseDelayMaxRatio);
    } else {
        delay =
            std::min(delay, decltype(buffer_)::kSize * (1.f - kSpeedModAmp));
    }
    // set new target speed
    targetSpeed_ = 1.f / delay * static_cast<float>(TapePosition::kUnity);

    setDrift(getDrift(), blockSize);

    maxBlockSizeWithDelay_ = std::min(
        maxBlockSize_, static_cast<int>(delay / (1.f - drift_ * kSpeedModAmp)));
}

void TapeDelay::setDrift(float drift, int blockSize)
{
    drift_ = drift;
    speedMod_.set({targetSpeed_ * drift * kSpeedModAmp},
                  1.f / static_cast<float>(blockSize));
}

void TapeDelay::setCutLowPass(float cutlowpass, int /*blockSize*/)
{
    cutlowpass_ = cutlowpass;
    auto freq   = cutlowpass_ * freqScale_;
    lpf_.setFreq({freq, freq});
}

void TapeDelay::setCutHiPass(float cuthighpass, int /*blockSize*/)
{
    cuthighpass_ = cuthighpass;
    auto freq    = cuthighpass_ * freqScale_;
    hpf_.setFreq({freq, freq});
}

void TapeDelay::setSaturation(float saturation, int blockSize)
{
    saturation_.set({saturation}, 1.f / static_cast<float>(blockSize));

    // update feedback
    setFeedback(getFeedback(), blockSize);
}
void TapeDelay::setFeedback(float feedback, int blockSize)
{
    feedback_ = feedback;

    // compensate saturation
    const auto saturation = saturation_.getTarget()[0];
    if (saturation >= 0.f) {
        feedback *= dsp::db2gain(-saturation / 2);
    } else // no increasing feedback if no saturation
    {
        feedback = std::min(feedback, 1.f);
    }

    feedbackCompensated_.set({feedback, feedback},
                             1.f / static_cast<float>(blockSize));
}
void TapeDelay::setDryWet(float drywet, int blockSize)
{
    drywet_            = drywet;
    float dry          = std::cos(dsp::constants<float>::pi_2 * drywet);
    float wet          = std::sin(dsp::constants<float>::pi_2 * drywet);
    float invBlockSize = 1.f / static_cast<float>(blockSize);
    dry_.set({dry, dry}, invBlockSize);
    wet_.set({wet, wet}, invBlockSize);
}

void TapeDelay::setMode(Mode mode, int blockSize)
{
    if (mode == kReverse) {
        tapReverse_.search(tapePos_);
    }
    switchTap(mode);
    setDelay(delay_, blockSize);
}

void TapeDelay::switchTap(Mode mode)
{
    oldMode_ = mode_;
    mode_    = mode;
    fadePos_ = kFadeSize - 1;
    tapId_ ^= 1;

    if (mode_ == kNormal) {
        tapTape_[tapId_].search(tapePos_);
    } else if (mode_ == kBackForth || mode == kReverse) {
        reverseDist_[tapId_] = 0;
        tapTape_[tapId_].reset(tapePos_);
    }

#ifdef TAPEDELAY_SWITCH_INDICATOR
    // set indicator
    switchIndicator_.store(true);
#endif
}

template <TapeDelay::Mode M, bool check = true, class Ctxt>
bool TapeDelay::read(Ctxt ctxt, int tapId,
                     TapeDelay::TapePosition::position_t speed)
{
    auto &tapTape = tapTape_[tapId];
    auto &x       = ctxt.getSignal();

    if constexpr (M == kNormal) {
        x = tapTape.read(ctxt, delayline_, tapePos_);
    } else if constexpr (M == kBackForth || M == kReverse) {
        auto &reverseDist = reverseDist_[tapId];

        reverseDist += speed + speed; // 2*speed

        // read tape
        x = tapTape.read<TapTape::kReverse>(ctxt, delayline_, tapePos_,
                                            reverseDist);

        if constexpr (M == kReverse) {
            auto xreverse = tapReverse_.read(ctxt, delaylineReverse_, tapePos_);
            inFor(x, k, i) { x[k][i] += xreverse[k][i]; }
        }

        // reach end of reverse
        if constexpr (check) {

            /* get reverse distance limit */
            constexpr auto kLimit = [] {
                if constexpr (M == kBackForth)
                    return static_cast<int>(TapePosition::kUnity);
                else if constexpr (M == kReverse)
                    return static_cast<int>(TapePosition::kUnity * 2);
            }();

            if (reverseDist > kLimit) {
                switchTap(M);
                reverseDist_[tapId_] =
                    reverseDist - kLimit + speed * kKernelSize * 2;
                return false;
            }
        }
    }
    return true;
}

TapeDelay::TapePosition::position_t TapeDelay::moveTape()
{
    // smooth speed;
    speed_ += (targetSpeed_ - speed_) * speedSmooth_;

    // speed modulation
    speedMod_.step();
    auto mod = speedLFO_.process()[0] * speedMod_.get()[0][0];

    // move tape
    auto speed = static_cast<TapePosition::position_t>(speed_ + mod);
    tapePos_.move(speed);

    return speed;
}

template <TapeDelay::Mode M, class Ctxt> int TapeDelay::readBlock(Ctxt ctxt)
{
    contextFor(ctxt)
    {
        auto speed = moveTape();

        if (!read<M>(c, tapId_, speed)) {
            return n + 1;
        }
    }
    return ctxt.getBlockSize();
}

void TapeDelay::process(const float *const *__restrict in,
                        float *const *__restrict out, int count)
{
    int blockSize = maxBlockSizeWithDelay_;

    const float *localin[] = {in[0], in[1]};
    float *localout[]      = {out[0], out[1]};

    while (count) {

        blockSize = std::min(blockSize, count);
        auto ctxt = dsp::BufferContext(x_, blockSize, buffer_);

        if (fadePos_ < 0) {
            // blockSize use in tape read part
            switch (mode_) {
            case kBackForth:
                blockSize = readBlock<kBackForth>(ctxt);
                break;
            case kReverse:
                blockSize = readBlock<kReverse>(ctxt);
                break;
            case kNormal:
            default:
                blockSize = readBlock<kNormal>(ctxt);
                break;
            }

            // update context
            ctxt.setBlockSize(blockSize);
        } else {
            blockSize = std::min(blockSize, fadePos_ + 1);
            ctxt.setBlockSize(blockSize);
            contextFor(ctxt)
            {
                auto speed = moveTape();

                auto &x = c.getSignal();
                switch (mode_) {
                case kBackForth:
                    read<kBackForth, false>(c, tapId_, speed);
                    break;
                case kReverse:
                    read<kReverse, false>(c, tapId_, speed);
                    break;
                case kNormal:
                default:
                    read<kNormal, false>(c, tapId_, speed);
                    break;
                }

                auto xIn = x;
                switch (oldMode_) {
                case kBackForth:
                    read<kBackForth, false>(c, tapId_ ^ 1, speed);
                    break;
                case kReverse:
                    read<kReverse, false>(c, tapId_ ^ 1, speed);
                    break;
                case kNormal:
                default:
                    read<kNormal, false>(c, tapId_ ^ 1, speed);
                    break;
                }

                auto fade = fadeLut[static_cast<size_t>(fadePos_)];
                inFor(x, k, i) { x[k][i] += fade * (xIn[k][i] - x[k][i]); }
                --fadePos_;
            }
        }

        // low pass filter
        contextFor(ctxt) { lpf_.process(c, lpfMem_); }

        // high pass filter
        contextFor(ctxt) { hpf_.process(c, hpfMem_); }

        // distortion
        if (saturation_.isActive() ||
            // check if its possible to vectorize given blocksize
            static_cast<size_t>(blockSize) % dsp::fSample<kN>::kVectorSize !=
                0) {
            float pregain  = 1.f;
            float postgain = 1.f;
            contextFor(ctxt)
            {
                saturation_.step();
                auto saturation = saturation_.get()[0][0];
                pregain         = dsp::db2gain(saturation);
                postgain =
                    dsp::db2gain(-saturation / (saturation > 0.f ? 2.f : 1.f));
                auto &x = c.getSignal();
                inFor(x, k, i)
                {
                    x[k][i] = postgain * dsp::tanh(x[k][i] * pregain);
                }
            }
            inFor(pregain_, k, i)
            {
                pregain_[k][i]  = pregain;
                postgain_[k][i] = postgain;
            }
        } else {
            contextFor(ctxt.vec())
            {
                auto &x = c.getSignal();
                inFor(x, k, i)
                {
                    x[k][i] *= pregain_[k][i];
                    x[k][i] = dsp::tanh(x[k][i]);
                    x[k][i] *= postgain_[k][i];
                }
            }
        }

        auto mode = mode_;
        contextFor(ctxt)
        {
            auto &loop = c.getSignal();
            decltype(c)::Type xin;

            inFor(xin, k, i) { xin[k][i] = *localin[i]++; }

            dry_.step();
            wet_.step();
            auto dry = dry_.get();
            auto wet = wet_.get();
            inFor(xin, k, i)
            {
                *localout[i]++ = xin[k][i] * dry[k][i] + loop[k][i] * wet[k][i];
            }

            feedbackCompensated_.step();
            auto feedback = feedbackCompensated_.get();
            decltype(c)::Type inloop;

            if (mode == kReverse) {
                delayline_.write(c, xin);
                inFor(xin, k, i) { inloop[k][i] = loop[k][i] * feedback[k][i]; }
                delaylineReverse_.write(c, inloop);
            } else {
                inFor(xin, k, i)
                {
                    inloop[k][i] = xin[k][i] + loop[k][i] * feedback[k][i];
                }
                delayline_.write(c, inloop);
            }
        }
        buffer_.setLimits();
        buffer_.nextBlock(ctxt);

        count -= blockSize;
    }

    dry_.reset();
    wet_.reset();
    feedbackCompensated_.reset();
    speedMod_.reset();
    if (saturation_.isActive()) {
        saturation_.reset();
        inFor(pregain_, k, i)
        {
            auto vecSize    = pregain_.size();
            pregain_[k][i]  = pregain_[vecSize - 1][i];
            postgain_[k][i] = postgain_[vecSize - 1][i];
        }
    }
}
} // namespace processors
