#include "TapeDelay.h"
#include "dsp/Utils.h"

namespace processors
{
TapeDelay::FadeLut TapeDelay::fadeLut;

void TapeDelay::update(float delay, float feedback, float cutlowpass,
                       float cuthighpass, float saturation, float drift,
                       Mode mode, float drywet)
{
    if (delay != delay_) {
        setDelay(delay);
    }
    if (drift != drift_) {
        setDrift(drift);
    }
    if (cutlowpass != cutlowpass_) {
        setCutLowPass(cutlowpass);
    }
    if (cuthighpass != cuthighpass_) {
        setCutHiPass(cuthighpass);
    }
    if (mode != mode_) {
        setMode(mode);
    }
    setSaturation(saturation);
    setFeedback(feedback);
    setDryWet(drywet);
}

void TapeDelay::setDelay(float delay)
{
    delay_ = delay;

    // limits;
    delay = std::max(delay, 0.f);
    if (mode_ == Reverse) {
        delay = std::min(delay,
                         MaxDelay * DefaultSampleRate * invSampleRate_ / 3.f);
    } else {
        delay = std::min(delay, MaxDelay * DefaultSampleRate * invSampleRate_);
    }
    // set new target speed
    targetSpeed_ = 1.f / delay_ * invSampleRate_ * TapePosition::Unity;

    setDrift(getDrift());
}

void TapeDelay::setDrift(float drift)
{
    drift_ = drift;
    speedMod_.set({targetSpeed_ * drift * speedModAmp}, invBlockSize_);
}

void TapeDelay::setCutLowPass(float cutlowpass)
{
    cutlowpass_ = cutlowpass;
    auto freq   = cutlowpass_ * freqScale_;
    lpf_.setFreq({freq, freq});
}

void TapeDelay::setCutHiPass(float cuthighpass)
{
    cuthighpass_ = cuthighpass;
    auto freq    = cuthighpass_ * freqScale_;
    hpf_.setFreq({freq, freq});
}

void TapeDelay::setSaturation(float saturation)
{
    saturation_.set({saturation}, invBlockSize_);

    // update feedback
    setFeedback(getFeedback());
}
void TapeDelay::setFeedback(float feedback)
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

    feedbackCompensated_.set({feedback, feedback}, invBlockSize_);
}
void TapeDelay::setDryWet(float drywet)
{
    drywet_.set({drywet, drywet}, invBlockSize_);
}

void TapeDelay::setMode(Mode mode)
{
    if (mode == Reverse) {
        tapReverse_.search(tapePos_);
    }
    switchTap(mode);
}

void TapeDelay::switchTap(Mode mode)
{
    oldMode_ = mode_;
    mode_    = mode;
    fadePos_ = FadeSize - 1;
    tapId_ ^= 1;

    if (mode_ == Normal) {
        tapTape_[tapId_].search(tapePos_);
    } else if (mode_ == BackForth || mode == Reverse) {
        reverseDist_[tapId_] = 0;
        tapTape_[tapId_].reset(tapePos_);
    }
}

template <TapeDelay::Mode M, bool check = true, class Ctxt>
bool TapeDelay::read(Ctxt ctxt, int tapId,
                     TapeDelay::TapePosition::position_t speed)
{
    auto &tapTape = tapTape_[tapId];
    auto &x       = ctxt.getSignal();

    if constexpr (M == Normal) {
        x = tapTape.read(ctxt, delayline_, tapePos_);
    } else if constexpr (M == BackForth || M == Reverse) {
        auto &reverseDist = reverseDist_[tapId];

        reverseDist += speed + speed; // 2*speed

        // read tape
        x = tapTape.read<TapTape::Reverse>(ctxt, delayline_, tapePos_,
                                           reverseDist);

        if constexpr (M == Reverse) {
            auto xreverse = tapReverse_.read(ctxt, delaylineReverse_, tapePos_);
            inFor(x, k, i) { x[k][i] += xreverse[k][i]; }
        }

        // reach end of reverse
        if constexpr (check) {

            /* get reverse distance limit */
            constexpr auto limit = [] {
                if constexpr (M == BackForth)
                    return static_cast<int>(TapePosition::Unity);
                else if constexpr (M == Reverse)
                    return static_cast<int>(TapePosition::Unity * 2);
            }();

            if (reverseDist > limit) {
                switchTap(M);
                reverseDist_[tapId_] =
                    reverseDist - limit + speed * KernelSize * 2;
                return false;
            }
        }
    }
    return true;
}

float TapeDelay::moveTape()
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
        float speed = moveTape();

        if (!read<M>(c, tapId_, speed)) {
            return n + 1;
        }
    }
    return ctxt.getBlockSize();
}

void TapeDelay::process(const float *const *__restrict in,
                        float *const *__restrict out, int count)
{
    int blockSizeOrig = std::min(count, MaxBlockSize);
    auto ctxt         = dsp::BufferContext(x_, blockSizeOrig, buffer_);

    const float *localin[] = {in[0], in[1]};
    float *localout[]      = {out[0], out[1]};

    while (count) {

        auto blockSize = std::min(blockSizeOrig, count);
        ctxt.setBlockSize(blockSize);

        if (fadePos_ < 0) {
            // blockSize use in tape read part
            switch (mode_) {
            case BackForth:
                blockSize = readBlock<BackForth>(ctxt);
                break;
            case Reverse:
                blockSize = readBlock<Reverse>(ctxt);
                break;
            case Normal:
            default:
                blockSize = readBlock<Normal>(ctxt);
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
                case BackForth:
                    read<BackForth, false>(c, tapId_, speed);
                    break;
                case Reverse:
                    read<Reverse, false>(c, tapId_, speed);
                    break;
                case Normal:
                default:
                    read<Normal, false>(c, tapId_, speed);
                    break;
                }

                auto xIn = x;
                switch (oldMode_) {
                case BackForth:
                    read<BackForth, false>(c, tapId_ ^ 1, speed);
                    break;
                case Reverse:
                    read<Reverse, false>(c, tapId_ ^ 1, speed);
                    break;
                case Normal:
                default:
                    read<Normal, false>(c, tapId_ ^ 1, speed);
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
            static_cast<size_t>(blockSize) % dsp::fSample<N>::VectorSize != 0) {
            float pregain  = 1.f;
            float postgain = 1.f;
            contextFor(ctxt)
            {
                saturation_.step();
                auto saturation = saturation_.get()[0][0];
                pregain         = dsp::db2gain(saturation);
                postgain = dsp::db2gain(-saturation / (saturation > 0 ? 2 : 1));
                auto &x  = c.getSignal();
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

            drywet_.step();
            auto drywet = drywet_.get();
            inFor(xin, k, i)
            {
                *localout[i]++ =
                    xin[k][i] + drywet[k][i] * (loop[k][i] - xin[k][i]);
            }

            feedbackCompensated_.step();
            auto feedback = feedbackCompensated_.get();
            decltype(c)::Type inloop;

            if (mode == Reverse) {
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

        ctxt.nextBlock();
        ctxt.save(buffer_);

        count -= blockSize;
    }

    drywet_.reset();
    feedbackCompensated_.reset();
    speedMod_.reset();
    if (saturation_.isActive()) {
        saturation_.reset();
        inFor(pregain_, k, i)
        {
            auto K          = pregain_.size();
            pregain_[k][i]  = pregain_[K - 1][i];
            postgain_[k][i] = postgain_[K - 1][i];
        }
    }
}
} // namespace processors
