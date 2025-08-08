#include "Adaptive_SC.h"
#include "dsp/AdaptiveFilter.h"
#include "dsp/AllPass.h"
#include "dsp/Context.h"

struct RLSUnit : public Unit {
    float fbufnum_;
    SndBuf *buf_;

    template <size_t Order, bool Warp> struct RLS {
        RLS(float lambda = 0.0995) : rls_{lambda} {}
        using aFilter = dsp::AdaptiveFilter<float, Order>;
        using DL = std::conditional_t<!Warp, dsp::CopyDelayLine<float, Order>,
                                      dsp::AllPassDelayLine<float, Order>>;
        dsp::RLS<float, Order> rls_;
        DL dlin;
    };
    void *rls_;
};

template <size_t Order, bool Warp>
static void rlsNext(RLSUnit *unit, int inNumSamples)
{
    float fbufnum = IN0(2);
    if (fbufnum != unit->fbufnum_) {
        uint32 bufnum = (int)fbufnum;
        World *world  = unit->mWorld;
        if (bufnum >= world->mNumSndBufs) bufnum = 0;
        unit->fbufnum_ = fbufnum;
        unit->buf_     = world->mSndBufs + bufnum;
    }

    const SndBuf *buf = unit->buf_;
    ACQUIRE_SNDBUF_SHARED(buf);

    auto *filter = (typename RLSUnit::RLS<Order, Warp>::aFilter *)buf->data;
    auto *rls    = (RLSUnit::RLS<Order, Warp> *)unit->rls_;

    if constexpr (Warp) {
        rls->dlin.setCoeff(IN0(4));
    }

    auto *in  = IN(0);
    auto *out = OUT(0);

    auto afilter = filter[inNumSamples];
    filter[0]    = afilter;

    for (int n = 0; n < inNumSamples; ++n) {
        auto x = *in;
        rls->rls_.process(dsp::Context(&x), rls->dlin, afilter);
        filter[n + 1] = afilter;
        *out          = x;
        ++in;
        ++out;
    }

    RELEASE_SNDBUF_SHARED(buf);
}

template <size_t Order, bool Warp> static void rlsAllocate(RLSUnit *unit)
{
    unit->rls_ = RTAlloc(unit->mWorld, sizeof(RLSUnit::RLS<Order, Warp>));
    new (unit->rls_) RLSUnit::RLS<Order, Warp>(IN0(3));

#define CMA ,
    SETCALC(rlsNext<Order CMA Warp>);
#undef CMA
}

template <bool Warp = false> static void rlsCtor(RLSUnit *unit)
{
    int order      = static_cast<int>(IN0(1));
    unit->fbufnum_ = -1e9f;
    switch (order) {
    case 2:
        rlsAllocate<2, Warp>(unit);
        break;
    case 3:
        rlsAllocate<3, Warp>(unit);
        break;
    case 4:
        rlsAllocate<4, Warp>(unit);
        break;
    case 5:
        rlsAllocate<5, Warp>(unit);
        break;
    case 6:
        rlsAllocate<6, Warp>(unit);
        break;
    case 7:
        rlsAllocate<7, Warp>(unit);
        break;
    case 8:
        rlsAllocate<8, Warp>(unit);
        break;
    default:
        return;
    }
}

static void rlsDtor(RLSUnit *unit) { RTFree(unit->mWorld, unit->rls_); }

struct ReconstructUnit : public Unit {
    float fbufnum_;
    SndBuf *buf_;

    template <size_t Order, bool Warp> struct Reconstruct;
    template <size_t Order> struct Reconstruct<Order, false> {
        using aFilter = dsp::AdaptiveFilter<float, Order>;
        using DL      = typename aFilter::DL;
        DL dlout;
    };
    template <size_t Order> struct Reconstruct<Order, true> {
        using aFilter = dsp::WarpedIIR<float, Order>;
        using State   = typename aFilter::State;
        State state;
    };
    void *reconstruct_;
};

template <size_t Order, bool Warp>
static void reconstructNext(ReconstructUnit *unit, int inNumSamples)
{
    using afilter = typename ReconstructUnit::Reconstruct<Order, Warp>::aFilter;

    float fbufnum = IN0(2);
    if (fbufnum != unit->fbufnum_) {
        uint32 bufnum = (int)fbufnum;
        World *world  = unit->mWorld;
        if (bufnum >= world->mNumSndBufs) bufnum = 0;
        unit->fbufnum_ = fbufnum;
        unit->buf_     = world->mSndBufs + bufnum;
    }

    const SndBuf *buf = unit->buf_;
    ACQUIRE_SNDBUF_SHARED(buf);

    auto *filter = (afilter *)buf->data;
    auto *reconstruct =
        (ReconstructUnit::Reconstruct<Order, Warp> *)unit->reconstruct_;

    for (int n = 0; n < inNumSamples; ++n) {
        auto x = IN(0)[n];
        if constexpr (!Warp)
            filter[n].reconstruct(dsp::Context(&x), reconstruct->dlout);
        else
            filter[n].reconstruct(dsp::Context(&x), reconstruct->state, IN0(3));
        OUT(0)[n] = x;
    }

    RELEASE_SNDBUF_SHARED(buf);
}

template <size_t Order, bool Warp>
static void reconstructAllocate(ReconstructUnit *unit)
{
    unit->reconstruct_ = RTAlloc(
        unit->mWorld, sizeof(ReconstructUnit::Reconstruct<Order, Warp>));
    new (unit->reconstruct_) ReconstructUnit::Reconstruct<Order, Warp>();

#define CMA ,
    SETCALC(reconstructNext<Order CMA Warp>);
#undef CMA
}

template <bool Warp = false> static void reconstructCtor(ReconstructUnit *unit)
{
    unit->fbufnum_ = -1e9f;
    int order      = static_cast<int>(IN0(1));
    switch (order) {
    case 2:
        reconstructAllocate<2, Warp>(unit);
        break;
    case 3:
        reconstructAllocate<3, Warp>(unit);
        break;
    case 4:
        reconstructAllocate<4, Warp>(unit);
        break;
    case 5:
        reconstructAllocate<5, Warp>(unit);
        break;
    case 6:
        reconstructAllocate<6, Warp>(unit);
        break;
    case 7:
        reconstructAllocate<7, Warp>(unit);
        break;
    case 8:
        reconstructAllocate<8, Warp>(unit);
        break;
    default:
        return;
    }
}

static void reconstructDtor(ReconstructUnit *unit)
{
    RTFree(unit->mWorld, unit->reconstruct_);
}

void loadAdaptive()
{
    (*ft->fDefineUnit)("RLS", sizeof(RLSUnit), (UnitCtorFunc)&rlsCtor<>,
                       (UnitDtorFunc)&rlsDtor, 0);
    (*ft->fDefineUnit)("RLSWarped", sizeof(RLSUnit),
                       (UnitCtorFunc)&rlsCtor<true>, (UnitDtorFunc)&rlsDtor, 0);
    (*ft->fDefineUnit)("AdaptiveReconstruct", sizeof(ReconstructUnit),
                       (UnitCtorFunc)reconstructCtor<>,
                       (UnitDtorFunc)&reconstructDtor, 0);
    (*ft->fDefineUnit)("AdaptiveReconstructWarped", sizeof(ReconstructUnit),
                       (UnitCtorFunc)reconstructCtor<true>,
                       (UnitDtorFunc)&reconstructDtor, 0);
}
