#include "TapeDelay.h"
#include "SC_Unit.h"

void TapeDelay::update(float delay) { speed_ = 1.f / delay / sampleRate_; }

void TapeDelay::process(float **__restrict in, float **__restrict out,
                        int count)
{
    auto *inl  = in[0];
    auto *inr  = in[1];
    auto *outl = out[0];
    auto *outr = out[1];

    int blockSize = std::min(count, MaxBlockSize);

    auto ctxt = dsp::BufferContext<dsp::Signal<N>, decltype(buffer_)>(
        nullptr, blockSize, buffer_);
    while (count) {

        contextFor(ctxt)
        {
            tapePos_.move(speed_);
            auto loop = tapTape_.read(c, delayline_, tapePos_);
            // auto loop = dsp::TapFix<24002>().read(c, delayline_);

            dsp::Signal<N>::Scalar x = {{{*(inl++), *(inr++)}}};

            inFor(x, k, i) { x[k][i] += loop[k][i] * 0.8f; }

            delayline_.write(c, x);

            *outl++ = x[0][0];
            *outr++ = x[0][1];
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
    float *in[2]  = {IN(1), IN(2)};
    float *out[2] = {OUT(0), OUT(1)};

    unit->tapedelay->update(IN0(0));
    unit->tapedelay->process(in, out, inNumSamples);
}

PluginLoad(SCTapeDelay)
{
    ft = inTable; // store pointer to InterfaceTable
    DefineDtorUnit(SCTapeDelay);
}
