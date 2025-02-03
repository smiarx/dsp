#include "Springs.h"
#include "SC_Unit.h"

static const Springs::MRs multirates{0.97f};

void Springs::update(float R, float freq, float Td, float T60)
{
    if (R != R_ || freq != freq_) {
        setPole(R, freq);
    }
    if(Td != Td_)
    {
        setTd(Td);
    }
    if(T60 != T60_)
    {
        setT60(T60);
    }
}

void Springs::setPole(float R, float freq)
{
        int M = sampleRate_ / 2.f / freq * 1.1f;
        M     = std::min(M, MaxDecimate);

        float f            = freq * freqScale_ * M;
        dsp::Signal<N> fs;
        dsp::Signal<N> Rs;
        for (int i = 0; i < N; ++i) {
            fs[i] = f*freqFactor[i];
            fs[i] = std::min(0.995f, fs[i]);
            fs[i] = std::max(0.005f, fs[i]);
            Rs[i] = R*freqFactor[i];
            Rs[i] = std::min(Rs[i], 1.f);
            Rs[i] = std::max(Rs[i], -1.f);
        }
        allpass_.setPole(Rs, fs);

        multirate_  = multirates.get(M);
        decimateId_ = 0;

        /* scipy.signal.ellip(10, 1, 60, 1.0, analog=True, output='sos') */
        const float analoglowpasssos[][2][3] = {
            {{1.00000000e-03, 0.00000000e+00, 1.43497622e-02},
             {1.00000000e+00, 4.79554974e-01, 1.41121071e-01}},
            {{1.00000000e+00, 0.00000000e+00, 2.24894555e+00},
             {1.00000000e+00, 2.72287841e-01, 5.26272885e-01}},
            {{1.00000000e+00, 0.00000000e+00, 1.33580569e+00},
             {1.00000000e+00, 1.11075112e-01, 8.24889759e-01}},
            {{1.00000000e+00, 0.00000000e+00, 1.12840768e+00},
             {1.00000000e+00, 3.84115770e-02, 9.56268315e-01}},
            {{1.00000000e+00, 0.00000000e+00, 1.07288063e+00},
             {1.00000000e+00, 9.05107247e-03, 9.99553018e-01}}};
        lowpass_.sosanalog(analoglowpasssos, fs);

        int oldM = M_;
        M_ = M;
        R_    = R;
        freq_ = freq;

        if(M != oldM)
        {
            setTd(Td_);
        }
}

void Springs::setTd(float Td)
{
    Td_ = Td;
    dsp::iSignal<N> loopEchoT;
    float sampleTd = Td *  sampleRate_ / M_;
    for(int i = 0; i < N; ++i)
    {
        loopTd_[i] = sampleTd * loopTdFactor[i];;
        loopModAmp_[i] = loopTd_[i]* loopModFactor[i];
        loopEchoT[i] = loopTd_[i] / 5.f;
        loopTd_[i] -= loopEchoT[i];
    }

    loopEcho_.setDelay(loopEchoT);

    setT60(T60_);
}

void Springs::setT60(float T60)
{
    T60_ = T60;
    loopGain_ = -powf(0.001, Td_/T60_);
}

void Springs::process(float **__restrict in, float **__restrict out, int count)
{
    auto *inl  = in[0];
    auto *inr  = in[1];
    auto *outl = out[0];
    auto *outr = out[1];

    int blockSize = std::min(count, MaxBlockSize);

    auto ctxt    = dsp::BufferContext(x_, blockSize, buffer_);
    auto ctxtdec = dsp::BufferContext(xdecimate_, blockSize, bufferDec_);
    while (count) {

        contextFor(ctxt)
        {
            auto &x = c.getIn();
            arrayFor(x, k)
            {
                auto mix = (*(inl++) + *(inr++)) * 0.5f;
                arrayFor(x[0], i) { x[k][i] = mix; }
            }
        }

        multirate_->decimate(ctxt, ctxtdec, dldecimate_, decimateId_);

        contextFor(ctxtdec)
        {
            auto &x   = c.getIn();

            LoopType looptap;
            auto mod = loopMod_.process();
            auto delay = loopTd_;
            for(int i = 0; i < N; ++i)
            {
                delay[i] += mod[i]*loopModAmp_[i];
            }
            looptap.setDelay(delay);

            auto loop = looptap.read(c, loopdl_);

            loopEchoDL_.write(c, loop);

            auto loopEchoVal = loopEcho_.read(c, loopEchoDL_);

            inFor(loop, k, i) {
                loopEchoVal[k][i] += loopEchoGain*(loop[k][i] - loopEchoVal[k][i]);
            }

            loopRippleDL_.write(c, loop);
            auto loopRippleVal = dsp::TapTail{}.read(c, loopRippleDL_);

            inFor(loop, k, i) {
                loopRippleVal[k][i] += loopRippleGain*(loopEchoVal[k][i] - loopRippleVal[k][i]);
            }


            inFor(x, k, i) {
                x[k][i] += loopRippleVal[k][i]*loopGain_;
            }
        }

        contextFor(ctxtdec)
        {
            for (int j = 0; j < CascadeL; ++j) {
                allpass_.process(c, allpassdl_[j]);
            }
        }

        contextFor(ctxtdec.vec())
        {
            auto &x = c.getIn();
            loopdl_.write(c, x);
        }

        contextFor(ctxtdec) { lowpass_.process(c, lowpassdl_); }

        decimateId_ =
            multirate_->interpolate(ctxtdec, ctxt, dlinterpolate_, decimateId_);

        contextFor(ctxt)
        {
            auto &x = c.getIn();
            {
                auto mix = 0.f;
                for (auto i = 0; i < N / 2; ++i) {
                    mix += x[0][i];
                }
                *(outl++) = mix / N * 2;
            }
            {
                auto mix = 0.f;
                for (auto i = N / 2; i < N; ++i) {
                    mix += x[0][i];
                }
                *(outr++) = mix / N * 2;
            }
        }

        ctxt.nextBlock();
        ctxtdec.nextBlock();

        count -= blockSize;
    }
    ctxt.save(buffer_);
    ctxtdec.save(bufferDec_);
}

#include "SC_PlugIn.h"

static InterfaceTable *ft;

struct SCSprings : public Unit {
    Springs *springs;
};

void SCSprings_Ctor(SCSprings *unit);
void SCSprings_Dtor(SCSprings *unit);
void SCSprings_next(SCSprings *unit, int inNumSamples);

void SCSprings_Ctor(SCSprings *unit)
{
    SETCALC(SCSprings_next);

    unit->springs = (Springs *)RTAlloc(unit->mWorld, sizeof(Springs));
    ClearUnitIfMemFailed(unit->springs);
    *unit->springs = Springs(SAMPLERATE);

    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
}

void SCSprings_Dtor(SCSprings *unit) { RTFree(unit->mWorld, unit->springs); }

void SCSprings_next(SCSprings *unit, int inNumSamples)
{
    float *in[2]  = {IN(4), IN(5)};
    float *out[2] = {OUT(0), OUT(1)};

    unit->springs->update(IN0(0), IN0(1), IN0(2), IN0(3));
    unit->springs->process(in, out, inNumSamples);
}

PluginLoad(SCSprings)
{
    ft = inTable; // store pointer to InterfaceTable
    DefineDtorUnit(SCSprings);
}
