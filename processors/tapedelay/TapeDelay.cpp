#include "TapeDelay.h"
#include "dsp/Utils.h"
#include "SC_Unit.h"

TapeDelay::FadeLut TapeDelay::fadeLut;

void TapeDelay::update(float delay, float feedback, float cutlowpass,
                       float cuthighpass, float saturation, float drift,
                       Mode mode, float drywet)
{
    if (delay != delay_ || drift_ != drift) {
        delay_ = delay;
        // set new target speed
        targetSpeed_ = 1.f / delay_ * invSampleRate_ * TapePosition::Unity;

        drift_ = drift;
        speedMod_.set({targetSpeed_ * drift * speedModAmp}, invBlockSize_);
    }
    if (cutlowpass != cutlowpass_) {
        cutlowpass_ = cutlowpass;
        auto freq   = cutlowpass_ * freqScale_;
        lpf_.setFreq({freq, freq});
    }
    if (cuthighpass != cuthighpass_) {
        cuthighpass_ = cuthighpass;
        auto freq    = cuthighpass_ * freqScale_;
        hpf_.setFreq({freq, freq});
    }
    if (mode != mode_) {
        switchTap(mode);
    }

    saturation_.set({saturation}, invBlockSize_);
    feedback_.set({feedback, feedback}, invBlockSize_);
    drywet_.set({drywet, drywet}, invBlockSize_);
}

void TapeDelay::switchTap(Mode mode)
{
    oldMode_ = mode_;
    mode_    = mode;
    fadePos_ = FadeSize - 1;
    tapId_ ^= 1;

    if (mode_ == Normal) {
    } else if (mode_ == BackForth) {
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
    } else if constexpr (M == BackForth) {
        auto &reverseDist = reverseDist_[tapId];
        reverseDist += 2 * speed;
        // reach end of reverse
        if constexpr (check) {
            constexpr auto limit = static_cast<int>(TapePosition::Unity);
            if (reverseDist > limit) {
                switchTap(BackForth);
                reverseDist_[tapId_] =
                    reverseDist - limit + speed * KernelSize * 2;
                return false;
            }
        }

        // read tape
        x = tapTape.read<TapTape::Reverse>(ctxt, delayline_, tapePos_,
                                           reverseDist);
    }
    return true;
}

template <TapeDelay::Mode M, class Ctxt> int TapeDelay::readBlock(Ctxt ctxt)
{
    contextFor(ctxt)
    {
        // smooth speed;
        speed_ += (targetSpeed_ - speed_) * speedSmooth_;

        // speed modulation
        speedMod_.step();
        auto mod = speedLFO_.process()[0] * speedMod_.get()[0][0];

        // move tape
        auto speed = static_cast<TapePosition::position_t>(speed_ + mod);
        tapePos_.move(speed);

        if (!read<M>(c, tapId_, speed)) {
            return n;
        }
    }
    return ctxt.getBlockSize();
}

void TapeDelay::process(float **__restrict in, float **__restrict out,
                        int count)
{
    int blockSizeOrig = std::min(count, MaxBlockSize);
    auto ctxt         = dsp::BufferContext(x_, blockSizeOrig, buffer_);

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
                // smooth speed;
                speed_ += (targetSpeed_ - speed_) * speedSmooth_;

                // speed modulation
                speedMod_.step();
                auto mod = speedLFO_.process()[0] * speedMod_.get()[0][0];

                // move tape
                auto speed =
                    static_cast<TapePosition::position_t>(speed_ + mod);
                tapePos_.move(speed);

                auto &x = c.getSignal();
                switch (mode_) {
                case BackForth:
                    read<BackForth, false>(c, tapId_, speed);
                    break;
                case Reverse:
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
                case Normal:
                default:
                    read<Normal, false>(c, tapId_ ^ 1, speed);
                    break;
                }

                auto fade = fadeLut[fadePos_];
                inFor(x, k, i) { x[k][i] += fade * (xIn[k][i] - x[k][i]); }
                --fadePos_;
            }
        }

        // low pass filter
        contextFor(ctxt) { lpf_.process(c, lpfMem_); }

        // high pass filter
        contextFor(ctxt) { hpf_.process(c, hpfMem_); }

        // distortion
        if (saturation_.isActive()) {
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

        contextFor(ctxt)
        {
            auto &loop = c.getSignal();
            decltype(c)::Type xin;

            inFor(xin, k, i) { xin[k][i] = *in[i]++; }

            drywet_.step();
            feedback_.step();
            auto drywet   = drywet_.get();
            auto feedback = feedback_.get();
            inFor(xin, k, i)
            {
                *out[i]++ = xin[k][i] + drywet[k][i] * (loop[k][i] - xin[k][i]);
            }

            dsp::fSample<N>::Scalar inloop;
            inFor(xin, k, i)
            {
                inloop[k][i] = xin[k][i] + loop[k][i] * feedback[k][i];
            }
            delayline_.write(c, inloop);
        }
        buffer_.setLimits();

        ctxt.nextBlock();
        ctxt.save(buffer_);

        count -= blockSize;
    }

    drywet_.reset();
    feedback_.reset();
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

#include "SC_PlugIn.h"

static InterfaceTable *ft;

struct SCTapeDelay : public Unit {
    TapeDelay *tapedelay;
};

void SCTapeDelay_Ctor(SCTapeDelay *unit);
void SCTapeDelay_Dtor(SCTapeDelay *unit);
void SCTapeDelay_next(SCTapeDelay *unit, int inNumSamples);

void SCTapeDelay_Ctor(SCTapeDelay *unit)
{
    SETCALC(SCTapeDelay_next);

    unit->tapedelay = (TapeDelay *)RTAlloc(unit->mWorld, sizeof(TapeDelay));
    ClearUnitIfMemFailed(unit->tapedelay);
    new (unit->tapedelay) TapeDelay(SAMPLERATE, BUFLENGTH);

    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
}

void SCTapeDelay_Dtor(SCTapeDelay *unit)
{
    RTFree(unit->mWorld, unit->tapedelay);
}

void SCTapeDelay_next(SCTapeDelay *unit, int inNumSamples)
{
    float *in[2]  = {IN(8), IN(9)};
    float *out[2] = {OUT(0), OUT(1)};

    unit->tapedelay->update(IN0(0), IN0(1), IN0(2), IN0(3), IN0(4), IN0(5),
                            static_cast<TapeDelay::Mode>(IN0(6)), IN0(7));
    unit->tapedelay->process(in, out, inNumSamples);
}

PluginLoad(SCTapeDelay)
{
    ft = inTable; // store pointer to InterfaceTable
    DefineDtorUnit(SCTapeDelay);
}
