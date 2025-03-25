#pragma once

#include "Delay.h"
#include "Signal.h"
#include <cmath>

namespace dsp
{
template <class... T> class List
{
};

template <size_t N, class Tap = TapTail> class AllPass
{
    /* Implement a Allpass filter with coeff a
     * template:
     *      N: vector size
     *      Tap: a Tap class for the recursive delayline
     *
     *  Allpass process needs a DelayLineType delayline to process
     *  internal state
     */
  public:
    constexpr AllPass() = default;
    template <class... Args>
    constexpr AllPass(const fData<N> &a, Args... args) : a_(a), tap_(args...)
    {
    }

    template <class Ctxt, class DL> void process(Ctxt c, DL &delayline) const
    {
        static_assert(Ctxt::VecSize <= DL::Length);

        auto &x = c.getSignal();
        auto sN = tap_.read(c, delayline);
        decltype(sN) s0;

        for (size_t k = 0; k < x.size(); ++k) {
            for (size_t i = 0; i < x[0].size(); ++i) {
                auto a = a_[i % N];

                s0[k][i] = x[k][i] - a * sN[k][i];
                x[k][i]  = a * s0[k][i] + sN[k][i];
            }
        }
        delayline.write(c, s0);
    }

    void setCoeff(fData<N> a) { a_ = a; }
    const auto &getCoeff() const { return a_; }
    void setDelay(fData<N> d) { tap_.setDelay(d); }
    void setDelay(iData<N> d) { tap_.setDelay(d); }

    const auto &getTap() const { return tap_; }

  private:
    fData<N> a_;
    Tap tap_;
};

template <size_t N, class TapOut = TapTail, class TapIn = TapTail>
class TapAllPass : public TapOut
{
  public:
    static constexpr auto minfdelay = 0.5f;

    constexpr TapAllPass() = default;
    template <class... Args>
    constexpr TapAllPass(const fData<N> &a, Args... args) :
        TapOut(args...), allpass_(a)
    {
    }

    void setDelay(fData<N> d)
    {
        fData<N> a;
        iData<N> id;
        for (int i = 0; i < N; ++i) {
            id[i]    = static_cast<int>(d[i] - minfdelay);
            float fd = d[i] - static_cast<float>(id[i]);
            a[i]     = (1 - fd) / (1 + fd);
        }
        allpass_.setCoeff(a);
        TapOut::setDelay(id);
    }

    template <class Ctxt, class DL> auto read(Ctxt c, DL &delayline) const
    {
        auto x = TapOut::read(c, delayline);
        c.setIn(x);
        allpass_.process(c, delayline.getInner());
        return x;
    }

  private:
    AllPass<N, TapIn> allpass_;
};

template <size_t N>
class TapAllPass<N, TapNoInterp<N>, TapTail> : TapNoInterp<N>
{
  public:
    static constexpr auto minfdelay = 0.618f;
    using TapOut                    = TapNoInterp<N>;
    using TapIn                     = TapTail;

    constexpr TapAllPass() = default;

    void setDelay(fData<N> d)
    {
        fData<N> a;
        iData<N> id;
        for (int i = 0; i < N; ++i) {
            id[i] = static_cast<int>(d[i] - minfdelay);

            float fd = d[i] - static_cast<float>(id[i]);
            // taylor approximation of (1-fd)/(1+fd)
            // -(fd-1)/2 + (fd-1)²/4 - (fd-1)³/8
            float z = (fd - 1.f) * 0.5f;
            a[i]    = z * (-1.f + z * (1.f - z));
        }
        allpass_.setCoeff(a);
        TapOut::setDelay(id);
    }

    template <class Ctxt, class DL> auto read(Ctxt c, DL &delayline) const
    {
        static_assert(Ctxt::VecSize == 1);
        typename Ctxt::Type x;
        typename Ctxt::Type x1;
        typename Ctxt::Type y;
        typename Ctxt::Type y1;

        for (int i = 0; i < N; i++) {
            auto id = TapOut::id_[i];
            assert(id + 1 <= DL::Length - Ctxt::VecSize + 1);

            auto &val = delayline.read(c, id);
            for (int k = 0; k < Ctxt::VecSize; ++k) x[k][i] = val[k][i];
            auto &val1 = delayline.read(c, id + 1);
            for (int k = 0; k < Ctxt::VecSize; ++k) x1[k][i] = val1[k][i];
        }

        y1 = allpass_.getTap().read(c, delayline.getInner());

        const auto &a = allpass_.getCoeff();
        for (int k = 0; k < Ctxt::VecSize; ++k) {
            for (int i = 0; i < N; i++) {
                y[k][i] = a[i] * (x[k][i] - y1[k][i]) + x1[k][i];
            }
        }
        delayline.getInner().write(c, y);
        return y;
    }

  private:
    AllPass<N, TapIn> allpass_;
};

template <size_t N> class AllPass2
{
    /* 2nd order allpass */

  public:
    struct State : std::array<fData<N>, 2> {
        State() : std::array<fData<N>, 2>{} {}
    };

    void setFreq(fData<N> freq)
    {
        for (size_t i = 0; i < N; ++i) {
            a_[1][i] = -std::cos(constants<float>::pi * freq[i]);
        }
    }

    void setRes(fData<N> R)
    {
        for (size_t i = 0; i < N; ++i) {
            assert(R[i] >= 0);
            a_[0][i] = (1.f - R[i]) / (1.f + R[i]);
        }
    }

    void setFreq(fData<N> freq, fData<N> R)
    {
        setFreq(freq);
        setRes(R);
    }

    void setCoeffs(fData<N> a0, fData<N> a1)
    {
        a_[0] = a0;
        a_[1] = a1;
    }

    template <class Ctxt, class State> void process(Ctxt c, State &state) const
    {
        auto &x = c.getSignal();
        auto s  = state;
        arrayFor(x, k)
        {
#pragma omp simd
            arrayFor(x[k], i)
            {
                auto &x0 = x[k][i];
                auto &s0 = s[0][i];
                auto &s1 = s[1][i];
                auto &a0 = a_[0][i];
                auto &a1 = a_[1][i];

                auto v0 = a0 * (x0 - s0);
                auto y0 = v0 + s0;
                auto x1 = v0 + x0;

                auto v1 = a1 * (x1 - s1);
                s0      = v1 + s1;
                s1      = v1 + x1;

                x[k][i] = y0;
            }
        }
        state = s;
    }

    __PROCESSBLOCK__

  private:
    fData<N> a_[2]; // allpass coeffs
};

// template <size_t N, bool EnergyPreserving = false> class AllPass2
//{
//     /* 2nd order allpass */
//
//   public:
//     struct State : std::array<fData<N>, 2> {
//         State() : std::array<fData<N>, 2>{} {}
//     };
//
//     void setFreq(fData<N> freq)
//     {
//         for (int i = 0; i < N; ++i) {
//             a_[1][i] = -cos(M_PIf * freq[i]);
//
//             if constexpr (EnergyPreserving) {
//                 //postgains_[1][i] = tan(M_PIf * freq[i] * 0.5f);
//                 postgains_[1][i] = sqrt((1.f+a_[1][i])/(1.f-a_[1][i]));
//                 //postgains_[1][i] = sqrt(1.f - a_[1][i]*a_[1][i]);
//             }
//         }
//         computeGains();
//     }
//
//     void setRes(fData<N> R)
//     {
//         for (int i = 0; i < N; ++i) {
//             assert(R[i] >= 0);
//             a_[0][i] = (1.f - R[i]) / (1.f + R[i]);
//
//             if constexpr (EnergyPreserving) {
//                 //postgains_[0][i] = sqrtf(R[i]);
//                 //postgains_[0][i] = sqrt((1.f+a_[0][i])/(1.f-a_[0][i]));
//                 postgains_[0][i] = sqrt(1.f - a_[0][i]*a_[0][i]);
//             }
//         }
//         computeGains();
//     }
//
//     void setFreq(fData<N> freq, fData<N> R)
//     {
//         setFreq(freq);
//         setRes(R);
//     }
//
//     void setCoeffs(fData<N> a0, fData<N> a1)
//     {
//         a_[0] = a0;
//         a_[1] = a1;
//     }
//
//     template <class Ctxt, class State> void process(Ctxt c, State &state)
//     //const
//     {
//         auto &x = c.getSignal();
//         auto s  = state;
//         arrayFor(x, k)
//         {
//             arrayFor(x[k], i)
//             {
//                 auto &x0 = x[k][i];
//                 auto &s0 = s[0][i];
//                 auto &s1 = s[1][i];
//                 auto &a0 = a_[0][i];
//                 auto &a1 = a_[1][i];
//
//                 auto g0 = sqrt((1.f-a0)/(1.f+a0));
//                 auto g1 = sqrt((1.f-a1)/(1.f+a1));
//                 //g0 = 1.f/g0;
//                 //g1 = 1.f/g1;
//
//                 x0 *= g0*g1;
//
//                 //auto v0 = a0 * (x0 - s0);
//                 auto v0 = a0 * (x0 - s0*g1/g1save);
//                 //auto y0 = v0 + s0;
//                 auto y0 = g1save/g1*v0 + s0;
//                 auto x1 = v0 + x0;
//
//                 //auto v0 = a0 * (x0 + s0);
//                 ////auto v0 = a0 * (x0 + s0*g1/g1save);
//                 //auto y0 = v0 + s0;
//                 ////auto y0 = g1save/g1*v0 + s0;
//                 //auto x1 = - v0 + x0;
//
//                 y0 /= g0*g1save;
//
//
//                 auto v1 = a1 * (x1 - s1);
//                 s0      = v1 + s1;
//                 s1      = v1 + x1;
//                 //auto v1 = a1 * (x1 + s1);
//                 //s0      = v1 + s1;
//                 //s1      = -v1 + x1;
//
//
//                 x[k][i] = y0;
//
//                 g1save = g1;
//
//
//             }
//         }
//         state = s;
//     }
//
//     __PROCESSBLOCK__;
//
//     template <class Ctxt> void pregain(Ctxt c) const
//     {
//         if constexpr (EnergyPreserving) {
//             auto &x = c.getSignal();
//             inFor(x, k, i) { x[k][i] *= pregain_[i]; }
//         }
//     }
//
//     template <class Ctxt> void postgain(Ctxt c) const
//     {
//         if constexpr (EnergyPreserving) {
//             auto &x = c.getSignal();
//             inFor(x, k, i) { x[k][i] *= postgain_[i]; }
//         }
//     }
//
//   private:
//     void computeGains()
//     {
//         // see
//         //
//         https://dafx2020.mdw.ac.at/proceedings/papers/DAFx2020_paper_59.pdf
//         // allow for energy preserving all pass
//
//         if constexpr (EnergyPreserving) {
//             for (size_t i = 0; i < N; ++i) {
//                 //postgain_[i] =
//                 sqrtf((1.f+a_[0][i])/(1.f-a_[0][i])*(1.f+a_[1][i])/(1.f-a_[1][i]));//postgains_[0][i]
//                 * postgains_[1][i]; postgain_[i] = postgains_[0][i] *
//                 postgains_[1][i]; pregain_[i]  = 1.f / postgain_[i];
//             }
//         }
//     }
//
//     // dummy empty element
//     struct Empty {
//     };
//
//     fData<N> a_[2]; // allpass coeffs
//     std::conditional_t<EnergyPreserving, fData<N>, Empty>
//         postgains_[2]; // post gain of each allpass
//     std::conditional_t<EnergyPreserving, fData<N>, Empty>
//         pregain_; // pre gain of nested all pass filter
//     std::conditional_t<EnergyPreserving, fData<N>, Empty>
//         postgain_; // post gain of nested all pass filter
//
//     float g1save{1.f};
// };
} // namespace dsp
