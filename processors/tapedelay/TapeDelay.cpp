#include "TapeDelay.h"
#include "dsp/Buffer.h"
#include "dsp/Context.h"
#include "dsp/FastMath.h"
#include "dsp/Utils.h"
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace processors
{
inline namespace DSP_ARCH_NAMESPACE
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
    decltype(ctxt.getInput()) x;

    bool continueBlock = true;

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
            x += xreverse;
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

                // write new tape read value
                ctxt.setOutput(x);
                continueBlock = false;
            }
        }
    }

    ctxt.setOutput(x);

    return continueBlock;
}

template <class Ctxt>
TapeDelay::TapePosition::position_t TapeDelay::moveTape(Ctxt ctxt)
{
    // smooth speed;
    speed_ += (targetSpeed_ - speed_) * speedSmooth_;

    // speed modulation
    auto mod = speedLFO_.process() * speedMod_.step(ctxt);

    // move tape
    auto speed = static_cast<TapePosition::position_t>(speed_ + mod);
    tapePos_.move(speed);

    return speed;
}

template <TapeDelay::Mode M, class Ctxt> int TapeDelay::readBlock(Ctxt ctxt)
{
    static_assert(!Ctxt::kUseVec);
    for (auto n = 0; n < ctxt.getBlockSize(); n += decltype(ctxt)::kIncrSize) {
        auto speed = moveTape(ctxt);

        if (!read<M>(ctxt, tapId_, speed)) {
            return n + 1;
        }
        ctxt.next();
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
            CTXTRUN(ctxt)
            {
                auto speed = moveTape(ctxt);

                switch (mode_) {
                case kBackForth:
                    read<kBackForth, false>(ctxt, tapId_, speed);
                    break;
                case kReverse:
                    read<kReverse, false>(ctxt, tapId_, speed);
                    break;
                case kNormal:
                default:
                    read<kNormal, false>(ctxt, tapId_, speed);
                    break;
                }
                auto xIn = ctxt.getInput();

                switch (oldMode_) {
                case kBackForth:
                    read<kBackForth, false>(ctxt, tapId_ ^ 1, speed);
                    break;
                case kReverse:
                    read<kReverse, false>(ctxt, tapId_ ^ 1, speed);
                    break;
                case kNormal:
                default:
                    read<kNormal, false>(ctxt, tapId_ ^ 1, speed);
                    break;
                }
                auto xOut = ctxt.getInput();

                auto fade = fadeLut[static_cast<size_t>(fadePos_)];
                auto x    = xOut + fade * (xIn - xOut);
                --fadePos_;

                ctxt.setOutput(x);
            };
        }

        CTXTRUN(ctxt)
        {
            // low pass filter
            lpf_.process(ctxt, lpfMem_);
            // high pass filter
            hpf_.process(ctxt, hpfMem_);

            // saturation
            if (saturation_.isActive()) {
                auto x      = ctxt.getInput();
                using ctxtT = decltype(x);

                auto saturation = saturation_.step(ctxt);
                auto pregain    = dsp::db2gain(saturation);
                auto postgain   = dsp::db2gain(
                    -ctxtT(saturation) /
                    dsp::blend(saturation > 0.f, ctxtT(2.), ctxtT(1.)));
                x = postgain_ * dsp::fasttanh(x * pregain_);
                ctxt.setOutput(x);

                pregain_  = pregain;
                postgain_ = postgain;
            } else {
                auto x = ctxt.getInput();
                x      = postgain_ * dsp::fasttanh(x * pregain_);
                ctxt.setOutput(x);
            }

            // write to delay buffers and output
            auto mode = mode_;
            {
                constexpr auto kIncrSize = decltype(ctxt)::kIncrSize;
                auto loop                = ctxt.getInput();
                mtype xin[kIncrSize];

                for (size_t j = 0; j < kIncrSize; ++j)
                    for (size_t i = 0; i < kN; ++i) xin[j][i] = *localin[i]++;

                auto dry = dry_.step(ctxt);
                auto wet = wet_.step(ctxt);

                mtype xout[kIncrSize];
                ctxt.store(xout[0], ctxt.load(xin[0]) * dry + loop * wet);
                for (size_t j = 0; j < kIncrSize; ++j)
                    for (size_t i = 0; i < kN; ++i) *localout[i]++ = xout[j][i];

                auto feedback = feedbackCompensated_.step(ctxt);

                if (mode == kReverse) {
                    delayline_.write(ctxt, ctxt.load(xin[0]));
                    auto inloop = loop * feedback;
                    delaylineReverse_.write(ctxt, inloop);
                } else {
                    auto inloop = ctxt.load(xin[0]) + loop * feedback;
                    delayline_.write(ctxt, inloop);
                }
            }
        };
        buffer_.nextBlock(ctxt);

        count -= blockSize;
    }

    dry_.reset();
    wet_.reset();
    feedbackCompensated_.reset();
    speedMod_.reset();
    saturation_.reset();
}

} // namespace DSP_ARCH_NAMESPACE
} // namespace processors
