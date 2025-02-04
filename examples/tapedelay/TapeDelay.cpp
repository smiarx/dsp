#include "TapeDelay.h"
#include "SC_Unit.h"

void TapeDelay::update(float delay, float feedback, float cutlowpass,
                       float cuthighpass, float drywet)
{
    if (delay != delay_) {
        delay_ = delay;
        // set new target speed
        targetSpeed_ = 1.f / delay_ * invSampleRate_ * TapePosition::Unity;
        speedMod_    = targetSpeed_ * 0.013f;
    }
    if (cutlowpass != cutlowpass_) {
        cutlowpass_ = cutlowpass;
        auto freq   = cutlowpass_ * freqScale_;
        lpf_.butterworthLP({freq, freq});
    }
    if (cuthighpass != cuthighpass_) {
        cuthighpass_ = cuthighpass;
        auto freq    = cuthighpass_ * freqScale_;
        printf("%f\n", freq);
        hpf_.butterworthHP({freq, freq});
    }

    feedback_ = feedback;
    drywet_   = drywet;
}

void TapeDelay::process(float **__restrict in, float **__restrict out,
                        int count)
{
    int blockSize = std::min(count, MaxBlockSize);

    auto ctxt = dsp::BufferContext(x_, blockSize, buffer_);
    while (count) {

        contextFor(ctxt)
        {
            auto &x = c.getIn();

            // smooth speed;
            speed_ += (targetSpeed_ - speed_) * speedSmooth_;

            // speed modulation
            auto mod = speedLFO_.process()[0] * speedMod_;

            // move tape
            tapePos_.move(static_cast<TapePosition::position_t>(speed_ + mod));
            x = tapTape_.read(c, delayline_, tapePos_);
        }

        // low pass filter
        contextFor(ctxt) { lpf_.process(c, lpfDL_); }
        // high pass filter
        contextFor(ctxt) { hpf_.process(c, hpfDL_); }

        contextFor(ctxt.vec())
        {
            auto &loop = c.getIn();
            decltype(c)::Type xin;

            inFor(xin, k, i) { xin[k][i] = *in[i]++; }
            inFor(xin, k, i)
            {
                *out[i]++ = xin[k][i] + drywet_ * (loop[k][i] - xin[k][i]);
            }

            inFor(xin, k, i) { xin[k][i] += loop[k][i] * feedback_; }

            delayline_.write(c, xin);
        }
        ctxt.nextBlock();
        ctxt.save(buffer_);

        count -= blockSize;
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
    new (unit->tapedelay) TapeDelay(SAMPLERATE);

    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
}

void SCTapeDelay_Dtor(SCTapeDelay *unit)
{
    RTFree(unit->mWorld, unit->tapedelay);
}

void SCTapeDelay_next(SCTapeDelay *unit, int inNumSamples)
{
    float *in[2]  = {IN(5), IN(6)};
    float *out[2] = {OUT(0), OUT(1)};

    unit->tapedelay->update(IN0(0), IN0(1), IN0(2), IN0(3), IN0(4));
    unit->tapedelay->process(in, out, inNumSamples);
}

PluginLoad(SCTapeDelay)
{
    ft = inTable; // store pointer to InterfaceTable
    DefineDtorUnit(SCTapeDelay);
}
