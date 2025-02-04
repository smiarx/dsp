#include "TapeDelay.h"
#include "SC_Unit.h"

void TapeDelay::update(float delay, float feedback, float drywet)
{
    if (delay != delay_) {
        delay_ = delay;
        // set new target speed
        targetSpeed_ = 1.f / delay_ * invSampleRate_ * TapePosition::Unity;
        speedMod_    = targetSpeed_ * 0.013f;
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
    float *in[2]  = {IN(3), IN(4)};
    float *out[2] = {OUT(0), OUT(1)};

    unit->tapedelay->update(IN0(0), IN0(1), IN0(2));
    unit->tapedelay->process(in, out, inNumSamples);
}

PluginLoad(SCTapeDelay)
{
    ft = inTable; // store pointer to InterfaceTable
    DefineDtorUnit(SCTapeDelay);
}
