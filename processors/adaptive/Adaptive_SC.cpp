#include "Adaptive_SC.h"
#include "dsp/AdaptiveFilter.h"
#include "dsp/AllPass.h"
#include "dsp/Context.h"

struct RLSUnit : public Unit {
    float fbufnum_;
    SndBuf *buf_;

    template <size_t Order, bool Warp> struct RLS {
        RLS(float lambda = 0.0995) : rls_{lambda} {}
        using aFilter =
            std::conditional_t<!Warp, dsp::AdaptiveFilter<float, Order>,
                               dsp::WarpedAdaptiveFilter<float, Order>>;
        dsp::RLS<float, Order> rls_;
        typename aFilter::AnalyzeState state_;
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

    auto *in  = IN(0);
    auto *out = OUT(0);

    auto afilter = filter[inNumSamples];
    if constexpr (Warp) {
        afilter.setWarpIn(IN0(4));
    }

    filter[0] = afilter;

    for (int n = 0; n < inNumSamples; ++n) {
        auto x = *in;
        rls->rls_.process(dsp::Context(&x), rls->state_, afilter);
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
    case 9:
        rlsAllocate<9, Warp>(unit);
        break;
    case 10:
        rlsAllocate<10, Warp>(unit);
        break;
    case 11:
        rlsAllocate<11, Warp>(unit);
        break;
    case 12:
        rlsAllocate<12, Warp>(unit);
        break;
    default:
        return;
    }
}

static void rlsDtor(RLSUnit *unit) { RTFree(unit->mWorld, unit->rls_); }

struct ReconstructUnit : public Unit {
    float fbufnum_;
    SndBuf *buf_;

    template <size_t Order, bool Warp> struct Reconstruct {
        using aFilter =
            std::conditional_t<!Warp, dsp::AdaptiveFilter<float, Order>,
                               dsp::WarpedAdaptiveFilter<float, Order>>;
        typename aFilter::ReconstructState state_;
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
        if constexpr (Warp) {
            filter[n].setWarpOut(IN0(3));
        }
        filter[n].reconstruct(dsp::Context(&x), reconstruct->state_);
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
    case 9:
        reconstructAllocate<9, Warp>(unit);
        break;
    case 10:
        reconstructAllocate<10, Warp>(unit);
        break;
    case 11:
        reconstructAllocate<11, Warp>(unit);
        break;
    case 12:
        reconstructAllocate<12, Warp>(unit);
        break;
    default:
        return;
    }
}

static void reconstructDtor(ReconstructUnit *unit)
{
    RTFree(unit->mWorld, unit->reconstruct_);
}

////////////////////////////////////////////////////////
struct FormantShiftUnit : public Unit {
    template <size_t Order, bool Warp> struct RLS {
        RLS(float lambda) : rls_(lambda) {}
        using aFilter =
            std::conditional_t<!Warp, dsp::AdaptiveFilter<float, Order>,
                               dsp::WarpedAdaptiveFilter<float, Order>>;
        aFilter filter_{};
        dsp::RLS<float, Order> rls_{};
        typename aFilter::AnalyzeState astate_;
        typename aFilter::ReconstructState rstate_;
    };
    void *rls_;
};

template <size_t Order, bool Warp>
static void formantShiftNext(FormantShiftUnit *unit, int inNumSamples)
{
    auto *rls = (FormantShiftUnit::RLS<Order, Warp> *)unit->rls_;

    if constexpr (Warp) {
        rls->filter_.setWarpIn(IN0(3));
        rls->filter_.setWarpIn(IN0(4));
    }

    auto *in  = IN(0);
    auto *out = OUT(0);

    rls->rls_.setForgetFactor(IN0(2));

    auto afilter = rls->filter_;
    for (int n = 0; n < inNumSamples; ++n) {
        auto x          = *in;
        auto origFilter = afilter;
        rls->rls_.process(dsp::Context(&x), rls->astate_, afilter);

        if constexpr (Warp) {
            afilter.compensateResidualAnalyze(dsp::Context(&x), rls->astate_);
            afilter.compensateResidualReconstruct(dsp::Context(&x),
                                                  rls->rstate_);
        }

        origFilter.reconstruct(dsp::Context(&x), rls->rstate_);

        *out = x;
        ++in;
        ++out;
    }
    rls->filter_ = afilter;
}

template <size_t Order, bool Warp>
static void formantShiftAllocate(FormantShiftUnit *unit)
{
    unit->rls_ =
        RTAlloc(unit->mWorld, sizeof(FormantShiftUnit::RLS<Order, Warp>));
    new (unit->rls_) FormantShiftUnit::RLS<Order, Warp>(IN0(2));

#define CMA ,
    SETCALC(formantShiftNext<Order CMA Warp>);
#undef CMA
}

template <bool Warp = false>
static void formantShiftCtor(FormantShiftUnit *unit)
{
    int order = static_cast<int>(IN0(1));
    switch (order) {
    case 2:
        formantShiftAllocate<2, Warp>(unit);
        break;
    case 3:
        formantShiftAllocate<3, Warp>(unit);
        break;
    case 4:
        formantShiftAllocate<4, Warp>(unit);
        break;
    case 5:
        formantShiftAllocate<5, Warp>(unit);
        break;
    case 6:
        formantShiftAllocate<6, Warp>(unit);
        break;
    case 7:
        formantShiftAllocate<7, Warp>(unit);
        break;
    case 8:
        formantShiftAllocate<8, Warp>(unit);
        break;
    case 9:
        formantShiftAllocate<9, Warp>(unit);
        break;
    case 10:
        formantShiftAllocate<10, Warp>(unit);
        break;
    default:
        return;
    }
}

static void formantShiftDtor(FormantShiftUnit *unit)
{
    RTFree(unit->mWorld, unit->rls_);
}

//==================================================
struct ANFUnit : public Unit {
    dsp::AdaptiveNotchFilter<float> anf_{};
    decltype(anf_)::State state_{};
};

static void anfNext(ANFUnit *unit, int inNumSamples)
{
    float *in  = IN(0);
    float *out = OUT(0);

    float rho = IN0(1);
    if (rho != unit->anf_.getRho()) unit->anf_.setRho(rho);

    float lambda = IN0(2);
    if (lambda != unit->anf_.getLambda()) unit->anf_.setLambda(lambda);

    for (int n = 0; n < inNumSamples; ++n) {
        auto x = in[n];
        unit->anf_.process(dsp::Context(&x), unit->state_);
        out[n] = unit->state_.getFreq() * static_cast<float>(SAMPLERATE) / 2.f;
    }
}

static void anfCtor(ANFUnit *unit)
{
    unit->state_ = {1.99, 1};
    SETCALC(anfNext);
}

//==================================================
void loadAdaptive()
{
    (*ft->fDefineUnit)("RLS", sizeof(RLSUnit), (UnitCtorFunc)&rlsCtor<>,
                       (UnitDtorFunc)&rlsDtor, 0);
    (*ft->fDefineUnit)("RLSWarped", sizeof(RLSUnit),
                       (UnitCtorFunc)&rlsCtor<true>, (UnitDtorFunc)&rlsDtor, 0);
    (*ft->fDefineUnit)("AdaptiveReconstruct", sizeof(ReconstructUnit),
                       (UnitCtorFunc)&reconstructCtor<>,
                       (UnitDtorFunc)&reconstructDtor, 0);
    (*ft->fDefineUnit)("AdaptiveReconstructWarped", sizeof(ReconstructUnit),
                       (UnitCtorFunc)&reconstructCtor<true>,
                       (UnitDtorFunc)&reconstructDtor, 0);
    (*ft->fDefineUnit)("FormantShift", sizeof(ReconstructUnit),
                       (UnitCtorFunc)&formantShiftCtor<true>,
                       (UnitDtorFunc)&formantShiftDtor, 0);
    (*ft->fDefineUnit)("AdaptiveNotchFilter", sizeof(ANFUnit),
                       (UnitCtorFunc)anfCtor, nullptr, 0);
}
