#include "SC_PlugIn.h"
#include "Reverb.h"

void Reverb::process(float **in, float **out, int count)
{
    assert(xblock != nullptr);

    int blockSize = std::min(count, MaxBlockSize);

    auto *inl  = in[0];
    auto *inr  = in[1];
    auto *outl = out[0];
    auto *outr = out[1];

    auto c1 = dsp::BufferContext(xblock, buffer1);
    c1.setBlockSize(blockSize);

    auto c2 = dsp::BufferContext(xblock, buffer2);
    c2.setBlockSize(blockSize);
    while (count) {

        processBlock(c1, {
            auto &x = c1.getIn();
            x[0]    = *(inl++) + *(inr++);
        });

        processBlock(c1, { allpass1.process(c1, apd1); });
        processBlock(c1, { allpass1.process(c1, apd2); });
        processBlock(c1, { allpass1.process(c1, apd3); });
        processBlock(c1, { allpass1.process(c1, apd4); });

        processBlock(c2, {
            auto &x = c2.getIn();
            x[1]    = x[0];
        });

        processBlock(c2, {
            constexpr dsp::TapFix<6261 COMMA 4460> tap;
            auto loop = tap.read(c2, feedbackd);
            loop = loop*fb;

            auto &x = c2.getIn();
            x[0] += loop[1];
            x[1] += loop[0];
        });

        processBlock(c2, {
                lpfilter.process(c2, lpfilterd);
                });

        processBlock(c2, { ap2n1.process(c2, dapd1); });
        processBlock(c2, { ap2n2.process(c2, dapd2); });
        processBlock(c2, { ap2n3.process(c2, dapd3); });
        processBlock(c2, { ap2n4.process(c2, dapd4); });

        processBlock(c2, {
                feedbackd.write(c2, c2.getIn());
        });

        processBlock(c2, {
                auto& x = c2.getIn();
               *(outl++) = x[0];
               *(outr++) = x[1];
        });

        c1.nextBlock();
        c2.nextBlock();
        count -= blockSize;
    }
    c1.save(buffer1);
    c2.save(buffer2);
}

static InterfaceTable *ft;

struct SCReverb : public Unit
{
    Reverb reverb;
};

void SCReverb_Ctor(SCReverb *unit);
void SCReverb_Dtor(SCReverb *unit);
void SCReverb_next(SCReverb *unit, int inNumSamples);

void SCReverb_Ctor(SCReverb *unit)
{
    SETCALC(SCReverb_next);

    unit->reverb = Reverb();

    auto* xblock =(dsp::Signal<2>*) RTAlloc(unit->mWorld, sizeof(dsp::Signal<2>)*Reverb::MaxBlockSize);
    ClearUnitIfMemFailed(xblock);
    unit->reverb.xblock = xblock;

    auto* buffer1 = (dsp::Signal<1>*) RTAlloc(unit->mWorld, sizeof(dsp::Signal<1>)*Reverb::Buffer1Size);
    ClearUnitIfMemFailed(buffer1);
    unit->reverb.buffer1.setBuffer(buffer1);

    auto* buffer2 = (dsp::Signal<2>*) RTAlloc(unit->mWorld, sizeof(dsp::Signal<2>)*Reverb::Buffer2Size);
    ClearUnitIfMemFailed(buffer2);
    unit->reverb.buffer2.setBuffer(buffer2);

    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
}

void SCReverb_Dtor(SCReverb *unit)
{
    RTFree(unit->mWorld,unit->reverb.xblock);
    RTFree(unit->mWorld,unit->reverb.buffer1.getBuffer());
    RTFree(unit->mWorld,unit->reverb.buffer2.getBuffer());
}

void SCReverb_next(SCReverb *unit, int inNumSamples)
{
    float* in[2] = {IN(0),IN(1)};
    float* out[2] = {OUT(0),OUT(1)};

    unit->reverb.process(in,out,inNumSamples);
}

PluginLoad(SCReverb)
{
    ft = inTable; // store pointer to InterfaceTable
    DefineDtorUnit(SCReverb);
}

