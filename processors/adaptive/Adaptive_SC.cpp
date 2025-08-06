#include "Adaptive_SC.h"
#include "dsp/AdaptiveFilter.h"
#include "dsp/Context.h"

struct RLSUnit : public Unit {
    float fbufnum_;
    SndBuf *buf_;

    template <size_t Order> struct RLS {
        using aFilter = dsp::AdaptiveFilter<float, Order>;
        using DL      = typename aFilter::DL;
        dsp::RLSDCD<float, Order> rls_{0.995};
        DL dlin;
    };
    void *rls_;
};

static void rlsCtor(RLSUnit *unit);

template <size_t Order> static void rlsNext(RLSUnit *unit, int inNumSamples)
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

    auto *filter = (typename RLSUnit::RLS<Order>::aFilter *)buf->data;
    auto *rls    = (RLSUnit::RLS<Order> *)unit->rls_;

    auto *in  = IN(0);
    auto *out = OUT(0);
    for (int n = 0; n < inNumSamples; ++n) {
        auto x = *in;
        rls->rls_.process(dsp::Context(&x), rls->dlin, *filter);
        *out = x;
        ++in;
        ++out;
    }

    RELEASE_SNDBUF_SHARED(buf);
}

template <size_t Order> static void rlsAllocate(RLSUnit *unit)
{
    unit->rls_ = RTAlloc(unit->mWorld, sizeof(RLSUnit::RLS<Order>));
    new (unit->rls_) RLSUnit::RLS<Order>();

    SETCALC(rlsNext<Order>);
}

void rlsCtor(RLSUnit *unit)
{
    int order      = static_cast<int>(IN0(1));
    unit->fbufnum_ = -1e9f;
    switch (order) {
    case 2:
        rlsAllocate<2>(unit);
        break;
    case 3:
        rlsAllocate<3>(unit);
        break;
    case 4:
        rlsAllocate<4>(unit);
        break;
    case 5:
        rlsAllocate<5>(unit);
        break;
    case 6:
        rlsAllocate<6>(unit);
        break;
    case 7:
        rlsAllocate<7>(unit);
        break;
    case 8:
        rlsAllocate<8>(unit);
        break;
    default:
        return;
    }
}

static void rlsDtor(RLSUnit *unit) { RTFree(unit->mWorld, unit->rls_); }

struct ReconstructUnit : public Unit {
    float fbufnum_;
    SndBuf *buf_;

    template <size_t Order> struct Reconstruct {
        using aFilter = dsp::AdaptiveFilter<float, Order>;
        using DL      = typename aFilter::DL;
        DL dlout;
    };
    void *reconstruct_;
};

static void reconstructCtor(ReconstructUnit *unit);

template <size_t Order>
static void reconstructNext(ReconstructUnit *unit, int inNumSamples)
{
    using afilter = dsp::AdaptiveFilter<float, Order>;

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
        (ReconstructUnit::Reconstruct<Order> *)unit->reconstruct_;

    for (int n = 0; n < inNumSamples; ++n) {
        auto x = IN(0)[n];
        filter->reconstruct(dsp::Context(&x), reconstruct->dlout);
        OUT(0)[n] = x;
    }

    RELEASE_SNDBUF_SHARED(buf);
}

template <size_t Order> static void reconstructAllocate(ReconstructUnit *unit)
{
    unit->reconstruct_ =
        RTAlloc(unit->mWorld, sizeof(ReconstructUnit::Reconstruct<Order>));
    new (unit->reconstruct_) ReconstructUnit::Reconstruct<Order>();

    SETCALC(reconstructNext<Order>);
}

static void reconstructCtor(ReconstructUnit *unit)
{
    unit->fbufnum_ = -1e9f;
    int order      = static_cast<int>(IN0(1));
    switch (order) {
    case 2:
        reconstructAllocate<2>(unit);
        break;
    case 3:
        reconstructAllocate<3>(unit);
        break;
    case 4:
        reconstructAllocate<4>(unit);
        break;
    case 5:
        reconstructAllocate<5>(unit);
        break;
    case 6:
        reconstructAllocate<6>(unit);
        break;
    case 7:
        reconstructAllocate<7>(unit);
        break;
    case 8:
        reconstructAllocate<8>(unit);
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
    (*ft->fDefineUnit)("RLS", sizeof(RLSUnit), (UnitCtorFunc)&rlsCtor,
                       (UnitDtorFunc)&rlsDtor, 0);
    (*ft->fDefineUnit)("AdaptiveReconstruct", sizeof(ReconstructUnit),
                       (UnitCtorFunc)reconstructCtor,
                       (UnitDtorFunc)&reconstructDtor, 0);
}
