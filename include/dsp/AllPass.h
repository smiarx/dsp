#pragma once

#include "Delay.h"
#include "Signal.h"
#include <cmath>

namespace dsp
{
template <class... T> class List
{
};

template <typename T, class Tap = TapTail> class AllPass
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
    constexpr AllPass(const T &a, Args... args) : a_(a), tap_(args...)
    {
    }

    template <class Ctxt, class DL> void process(Ctxt c, DL &delayline) const
    {
        static_assert(Ctxt::kIncrSize <= DL::kLength);

        auto x  = c.getInput();
        auto sN = tap_.read(c, delayline);

        auto s0 = x - a_ * sN;
        x       = a_ * s0 + sN;

        delayline.write(c, s0);
        c.setOutput(x);
    }

    void setCoeff(const T &a) { a_ = a; }
    [[nodiscard]] const auto &getCoeff() const { return a_; }
    void setDelay(const T &d) { tap_.setDelay(d); }
    void setDelay(const intType<T> &d) { tap_.setDelay(d); }

    [[nodiscard]] const auto &getTap() const { return tap_; }

  private:
    T a_;
    Tap tap_;
};

template <typename T, class TapOut = TapTail, class TapIn = TapTail>
class TapAllPass : public TapOut
{
  public:
    static constexpr auto kMinfdelay = 0.5f;

    constexpr TapAllPass() = default;
    template <class... Args>
    constexpr TapAllPass(const T &a, Args... args) :
        TapOut(args...), allpass_(a)
    {
    }

    void setDelay(const T &d)
    {
        auto id = toInt(load(d)) - kMinfdelay;
        auto fd = d - id;
        auto a  = (1 - fd) / (1 + fd);

        allpass_.setCoeff(a);
        TapOut::setDelay(id);
    }

    template <class Ctxt, class DL> auto read(Ctxt c, DL &delayline) const
    {
        T x = TapOut::read(c, delayline);

        c.setData(&x);
        allpass_.process(c, delayline.getInner());

        return c.getInput();
    }

  private:
    AllPass<T, TapIn> allpass_;
};

template <typename T>
class TapAllPass<T, TapNoInterp<T>, TapTail> : TapNoInterp<T>
{
  public:
    static constexpr auto kMinfdelay = 0.618f;
    using TapOut                     = TapNoInterp<T>;
    using TapIn                      = TapTail;

    constexpr TapAllPass() = default;

    void setDelay(const T &d)
    {
        auto id = toInt(load(d) - kMinfdelay);
        auto fd = d - id;
        // taylor approximation of (1-fd)/(1+fd)
        // -(fd-1)/2 + (fd-1)²/4 - (fd-1)³/8
        auto z = (fd - T(1)) * T(0.5);
        auto a = z * (T(-1) + z * (T(1) - z));
        allpass_.setCoeff(a);
        TapOut::setDelay(id);
    }

    template <class Ctxt, class DL> auto read(Ctxt c, DL &delayline) const
    {
        static_assert(Ctxt::IncrSize == 1);

        auto id = load(TapOut::id_);
        auto a  = allpass_.getCoeff();

        assert(all(load(id) + 1 <= DL::Length - Ctxt::IncrSize + 1));

        auto x  = TapNoInterp<T>{id}.read(c, delayline);
        auto x1 = TapNoInterp<T>{id + 1}.read(c, delayline);

        auto y1 = allpass_.getTap().read(c, delayline.getInner());
        auto y  = a * (x - y1) + x1;

        delayline.getInner().write(c, y);
        return y;
    }

  private:
    AllPass<T, TapIn> allpass_;
};

template <typename T> class AllPass2
{
    /* 2nd order allpass */

  public:
    struct State : std::array<T, 2> {
        State() : std::array<T, 2>{} {}
    };

    void setFreq(const T &freq) { a_[1] = -cos(constants<T>::pi * load(freq)); }

    void setRes(const T &r)
    {
        auto vr = load(r);
        assert(all(vr >= 0));
        a_[0] = (-vr + 1) / (vr + 1);
    }

    void setFreq(const T &freq, const T &r)
    {
        setFreq(freq);
        setRes(r);
    }

    void setCoeffs(const T &a0, const T &a1)
    {
        a_[0] = a0;
        a_[1] = a1;
    }

    template <class Ctxt, class State> void process(Ctxt c, State &state) const
    {
        static_assert(Ctxt::kIncrSize == 1, "Cannot vectorize AllPass2");
        auto &s = state;
        auto x0 = c.getInput();

        auto v0 = a_[0] * (x0 - s[0]);
        auto y0 = v0 + s[0];
        auto x1 = v0 + x0;

        auto v1 = a_[1] * (x1 - s[1]);
        s[0]    = v1 + s[1];
        s[1]    = v1 + x1;

        c.setOutput(y0);
    }

    PROCESSBLOCK_

  private:
    T a_[2]; // allpass coeffs
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
