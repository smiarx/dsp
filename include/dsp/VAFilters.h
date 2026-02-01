#pragma once

#include "Context.h"
#include "FastMath.h"

namespace dsp::va
{

/* Various implementation from
 * https://archive.org/details/the-art-of-va-filter-design-rev.-2.1.2
 */

enum FilterType {
    kLowPass,
    kHighPass,
    kAllPass,
    kBandPass,
    kNotch,
};

template <typename T> static constexpr auto warpGain(const T &freq)
{
    return tan(dsp::constants<T>::pi * load(freq) * 0.5);
}

template <typename T, FilterType FT = kLowPass> class OnePole
{
  public:
    struct State : public std::array<T, 1> {
        State() : std::array<T, 1>{} {}
    };

    OnePole(const T &freq)
    {
        static_assert(FT == kLowPass || FT == kHighPass || FT == kAllPass);
        setFreq(freq);
    }
    OnePole() = default;

    void setFreq(const T &freq)
    {
        auto gain = warpGain(freq);
        gain_     = gain / (load(T(1)) + gain);
    }

    template <class Ctxt> void process(Ctxt c, State &state) const
    {
        auto x = c.getInput();
        auto s = load(state[0]);

        auto v  = (x - s) * gain_;
        auto lp = v + s;
        s       = lp + v;

        if constexpr (FT == kLowPass) {
            x = lp;
        } else {
            auto hp = x - lp;
            if constexpr (FT == kHighPass) {
                x = hp;
            } else if constexpr (FT == kAllPass) {
                auto ap = lp - hp;
                x       = ap;
            }
        }

        c.setOutput(x);
        state[0] = s;
    }

    PROCESSBLOCK_

    [[nodiscard]] auto getGain() const
    {
        auto gain = load(gain_);
        if constexpr (FT == kLowPass) {
            return gain;
        } else if constexpr (FT == kHighPass) {
            gain = T(1) - gain;
        } else if constexpr (FT == kAllPass) {
            gain += gain;
        }
        return gain;
    }

  private:
    T gain_;
};

template <typename T, FilterType FT = kLowPass> class SVF
{
  public:
    static constexpr bool kNormalizedBandPass =
        FT == kBandPass || FT == kAllPass || FT == kNotch;

    struct State : public std::array<T, 2> {
        State() : std::array<T, 2>{} {}
    };

    SVF(const T &freq, const T &res) { setFreq(freq, res); }
    SVF() = default;

    void setFreq(const T &freq, const T &res)
    {
        gain_ = warpGain(freq);
        setRes(res);
    }
    void setFreq(const T &freq)
    {
        static constexpr auto kDefaultRes = dsp::constants<T>::sqrt1_2;
        T res                             = kDefaultRes;
        setFreq(freq, res);
    }

    void setRes(const T &res)
    {
        auto r    = load(res);
        auto gain = load(gain_);
        assert(all(r >= 0));

        denominator_ = T(1) / (T(1) + gain * (T(2) * r + gain));

        if constexpr (FT == kHighPass) {
            gains1_ = r * 2 + gain;
        }
        if constexpr (kNormalizedBandPass) {
            inputGain_ = r * 2;
            if constexpr (FT == kAllPass) {
                inputGain_ = load(inputGain_) * 2;
            }
        }
    }

    /* -3db bandwidth in octave */
    template <bool PreWarp = false> void setBandWidth(const T &bandwith)
    {
        static_assert(FT == kBandPass);
        using bT = baseType<T>;

        auto gain = load(gain_);
        auto bw   = load(bandwith);

        if constexpr (PreWarp) {
            auto freqbwm = gain * pow(bT(2), -bw * 0.5);
            auto freqbwp = gain * pow(bT(2), bw * 0.5);
            bw           = log(warpGain(freqbwp) / warpGain(freqbwm)) / log(2.);
        }

        auto res = pow(bT(2), bw * 0.5) - pow(bT(2), -bw * 0.5) * 0.5;
        setRes(res);
    }

    template <class Ctxt> void process(Ctxt c, State &state) const
    {
        static_assert(!Ctxt::kUseVec);

        auto x  = c.getInput();
        auto s1 = load(state[0]);
        auto s2 = load(state[1]);
        auto y  = x;

        if constexpr (kNormalizedBandPass) {
            // normalize
            x *= inputGain_;
        }

        if constexpr (FT == kHighPass) {
            auto hp = (x - gains1_ * s1 - s2) * denominator_;
            auto v1 = gain_ * hp;
            auto bp = v1 + s1;
            s1      = bp + v1; // first integrator
            auto v2 = gain_ * bp;
            auto lp = v2 + s2;
            s2      = lp + v2; // second integrator
            y       = hp;
        } else {
            auto bp = (gain_ * (x - s2) + s1) * denominator_;
            if constexpr (kNormalizedBandPass) {
                auto bp2 = bp + bp; // first integrator
                s1       = bp2 - s1;
                auto v22 = gain_ * bp2; // second integrator
                s2       = v22 + s2;
                if constexpr (FT == kBandPass) {
                    y = bp;
                } else if constexpr (FT == kAllPass || FT == kNotch) {
                    y = y - bp;
                }
            } else if (FT == kLowPass) {
                auto v1 = bp - s1; // first integrator
                s1      = bp + v1;
                auto v2 = gain_ * bp; // secondintegrator
                auto lp = v2 + s2;
                s2      = lp + v2;
                y       = lp;
            }
        }

        state[0] = s1;
        state[1] = s2;
        c.setOutput(y);
    }

    PROCESSBLOCK_

  private:
    T gain_;
    T denominator_;

    class Empty
    {
    };
    std::conditional_t<FT == kHighPass, T, Empty> gains1_;
    std::conditional_t<kNormalizedBandPass, T, Empty> inputGain_;
};

template <typename T, FilterType FT = kLowPass> class Ladder
{
  public:
    using OP    = OnePole<T, FT>;
    using State = std::array<typename OP::State, 4>;

    Ladder() = default;
    Ladder(const T &freq, const T &res) : onepole_(freq) { setFreq(freq, res); }

    void setFreq(const T &freq, const T &res = 0.f)
    {
        onepole_.setFreq(freq);
        setRes(res);
    }
    void setRes(const T &res)
    {
        auto gain    = onepole_.getGain();
        auto gainSq  = gain * gain;
        res_         = res;
        denominator_ = T(1) / (T(1) + gainSq * gainSq * res);
    }

    template <class Ctxt> void process(Ctxt c, State &state) const
    {
        auto x = c.getInput();

        auto g = onepole_.getGain();

        auto s1 = load(state[0][0]);
        auto s2 = load(state[1][0]);
        auto s3 = load(state[2][0]);
        auto s4 = load(state[3][0]);
        auto S  = s4 + g * (s3 + g * (s2 + g * s1));
        // https://www.kvraudio.com/forum/viewtopic.php?f=33&t=571909
        if constexpr (FT == kLowPass) {
            S *= (T(1) - g);
        } else if constexpr (FT == kHighPass) {
            S *= -g;
        }

        x = (x - res_ * S) * denominator_;

        c.setOutput(x);
        for (size_t j = 0; j < 4; ++j) {
            onepole_.process(c, state[j]);
        }
    }

    PROCESSBLOCK_

  private:
    T res_;
    T denominator_;
    OP onepole_;
};

template <typename T>
struct SparseJacobian
{
    T a,b,c,d,ap,bp,cp,dp;
};

template <typename T>
static inline auto solveNewton4Stages(SparseJacobian<T>& j, T*x, T *fx)
{
    auto a = j.a;
    auto b = j.b;
    auto c = j.c;
    auto d = j.d;
    auto ap = j.ap;
    auto bp = j.bp;
    auto cp = j.cp;
    auto dp = j.dp;

    T xmfx[4], solve[4];
    for (int i = 0; i < 4; ++i)
        xmfx[i] = x[i] - fx[i];

    T* v = xmfx;

    auto div = a * b * c * d - ap * bp * cp * dp;
    solve[0] = b * (c * (d * v[0] + dp * v[3]) + cp * dp * v[2]) +
               bp * cp * dp * v[1];
    solve[1] = c * (d * (a * v[1] + ap * v[0]) + dp * ap * v[3]) +
               cp * ap * dp * v[2];
    solve[2] = d * (a * (b * v[2] + bp * v[1]) + ap * bp * v[0]) +
               dp * ap * bp * v[3];
    solve[3] = a * (b * (c * v[3] + cp * v[2]) + bp * cp * v[1]) +
               ap * bp * cp * v[4];

    for (int i = 0; i < 4; ++i)
        x[i] -= solve[i] / div;
}

template <typename T, FilterType FT = kLowPass> class LadderOTA
{
  public:
    using State = struct State {
        T s_[4]{};
        T u_[4]{};
        T tanhu_[4]{};
    };

    LadderOTA() = default;
    LadderOTA(const T &freq, const T &res) { setFreq(freq, res); }

    void setFreq(const T &freq, const T &res = 0.f)
    {
        gain_ = warpGain(freq);
        setRes(res);
    }
    void setRes(const T &res) { res_ = res; }

    template <class Ctxt> void solveHP(Ctxt ctxt, State &state) const
    {
        auto x = ctxt.getInput();

        T R[4];
        for (int i = 0; i < 4; ++i)
            R[i] = 1 + gain_ * (1 - state.tanhu_[i] * state.tanhu_[i]);

        auto a = R[0], b = R[1], c = R[2], d = R[3];
        auto e = -res_;

        T gtanhus[4];
        T fu[4];
        for (int i = 0; i < 4; ++i) {
            gtanhus[i] = gain_ * state.tanhu_[i] + state.s_[i];
            fu[i]      = +state.u_[i] + gtanhus[i];
        }
        fu[0] -= x - res_ * state.u_[3];
        fu[1] -= state.u_[0];
        fu[2] -= state.u_[1];
        fu[3] -= state.u_[2];

        T efu[4];
        for (int i = 0; i < 4; ++i) efu[i] = e * fu[i];

        T solve[4];
        solve[0] = efu[1] + b * (efu[2] + c * (efu[3] + d * fu[0]));
        solve[1] = efu[2] + c * (efu[3] + d * (fu[0] + a * fu[1]));
        solve[2] = efu[3] + d * (fu[0] + a * (fu[1] + b * fu[2]));
        solve[3] = fu[0] + a * (fu[1] + b * (fu[2] + c * fu[3]));

        auto div = a * b * c * d - e;
        for (int i = 0; i < 4; ++i) {
            state.u_[i] -= solve[i] / div;
            state.tanhu_[i] = dsp::tanh(state.u_[i]);
        }
    }

    template <class Ctxt> void solveLP(Ctxt ctxt, State &state) const
    {
        auto x = ctxt.getInput();

        T gtanhu[4];
        for (int i = 0; i < 4; ++i) {
            gtanhu[i] = gain_ * state.tanhu_[i];
        }

        // derivation of // g*tanh(u) + s
        T deriv[4];
        for (int i = 0; i < 4; ++i)
            deriv[i] = gain_ - gtanhu[i] * state.tanhu_[i];

        // jacobian matrix coefficient
        auto a = 1 + deriv[0], b = 1 + deriv[1], c = 1 + deriv[2],
             d  = 1 + deriv[3];
        auto ap = deriv[0], bp = deriv[1], cp = deriv[2], dp = -res_ * deriv[3];

        // compute the value of u - f(u)
        T gtanhus[4];
        for (int i = 0; i < 4; ++i) {
            gtanhus[i] = gtanhu[i] + state.s_[i];
        }
        T fu[4];
        for (int i = 0; i < 4; ++i) {
            auto input = i == 0 ? x - res_ * gtanhus[3] : gtanhus[i - 1];
            fu[i]      = input - gtanhus[i];
        }

        // solve the newton value by computing direct invert matrix
        // multiplication
        // should be able to use fewer multiplications
        SparseJacobian<T> jacobian = {a,b,c,d,ap,bp,cp,dp};
        solveNewton4Stages(jacobian,state.u_,fu);

        for (int i = 0; i < 4; ++i) {
            state.tanhu_[i] = dsp::fasttanh(state.u_[i]);
        }
    }

    // template <class Ctxt> void solveLP2(Ctxt ctxt, State &state) const
    //{
    //     auto x = ctxt.getInput();
    //     auto& y = state.u_;

    //    T tanhy[4];
    //    for (int i = 0; i < 4; ++i) {
    //        auto input = i == 0 ? x - res_*y[3] : y[i-1];
    //        tanhy[i] = dsp::fasttanh(input-y[i]);
    //    }

    //    // derivation of // g*tanh(u) + s
    //    T deriv[4];
    //    for (int i = 0; i < 4; ++i)
    //        deriv[i] = gain_*(1 - tanhy[i]*tanhy[i]);

    //    // jacobian matrix coefficient
    //    auto a = 1 + deriv[0], b = 1 + deriv[1], c = 1 + deriv[2],
    //         d  = 1 + deriv[3];
    //    auto ap = deriv[0], bp = deriv[1], cp = deriv[2], dp = -res_ *
    //    deriv[3];

    //    // compute the value of u - f(u)
    //    T gtanhys[4];
    //    for (int i = 0; i < 4; ++i) {
    //        gtanhys[i] = gain_*tanhy[i] + state.s_[i];
    //    }
    //    T fy[4];
    //    for (int i = 0; i < 4; ++i) {
    //        fy[i]      = y[i] - gtanhys[i];
    //    }

    //    // solve the newton value by computing direct invert matrix
    //    // multiplication
    //    // should be able to use fewer multiplications
    //    T solve[4];
    //    solve[0] = b * (c * (d * fy[0] + dp * fy[3]) + cp * dp * fy[2]) +
    //               bp * cp * dp * fy[1];
    //    solve[1] = c * (d * (a * fy[1] + ap * fy[0]) + dp * ap * fy[3]) +
    //               cp * ap * dp * fy[2];
    //    solve[2] = d * (a * (b * fy[2] + bp * fy[1]) + ap * bp * fy[0]) +
    //               dp * ap * bp * fy[3];
    //    solve[3] = a * (b * (c * fy[3] + cp * fy[2]) + bp * cp * fy[1]) +
    //               ap * bp * cp * fy[4];

    //    auto div = a * b * c * d - ap * bp * cp * dp;
    //    for (int i = 0; i < 4; ++i) {
    //        y[i] -= solve[i] / div;
    //    }
    //}

    template <class Ctxt> void process(Ctxt ctxt, State &state) const
    {
        // if constexpr (FT == kLowPass)
        //{
        //     auto& y = state.u_;
        //     solveLP2(ctxt, state);
        //     for (int i = 0; i < 4; ++i) {
        //         state.s_[i] = y[i]+y[i] - state.s_[i];
        //     }
        //     ctxt.setOutput(y[3]);
        // }

        if constexpr (FT == kLowPass) {
            solveLP(ctxt, state);
        } else if constexpr (FT == kHighPass)
            solveHP(ctxt, state);

        T y;
        for (int i = 0; i < 4; ++i) {
            auto v      = gain_ * state.tanhu_[i];
            y           = v + state.s_[i];
            state.s_[i] = y + v;
        }

        if constexpr (FT == kLowPass) ctxt.setOutput(y);
        else
            ctxt.setOutput(state.u_[3]);
    }

    PROCESSBLOCK_

  private:
    T gain_;
    T res_;
};

template <typename T, FilterType FT = kLowPass> class LadderTransistor
{
  public:
    using State = struct State {
        T s_[4]{};
        T u_[4]{};
    };

    LadderTransistor() = default;
    LadderTransistor(const T &freq, const T &res) { setFreq(freq, res); }

    void setFreq(const T &freq, const T &res = 0.f)
    {
        gain_ = warpGain(freq);
        setRes(res);
    }
    void setRes(const T &res) { res_ = res; }

    template <class Ctxt> void solveLP(Ctxt ctxt, State &state) const
    {
        auto x  = ctxt.getInput();
        auto &u = state.u_;
        auto &s = state.s_;

        T tanhus[4];
        for (int i = 0; i < 4; ++i) {
            tanhus[i] = dsp::fasttanh(gain_ * u[i] + s[i]);
        }
        T tanhinput;
        if constexpr (FT == kLowPass)
            tanhinput = dsp::fasttanh(x - res_ * (gain_ * u[3] + s[3]));
        else if (FT == kHighPass)
            tanhinput = dsp::fasttanh(x - res_ * u[3]);

        // derivation of tanh(g*u+s)
        T deriv[4];
        for (int i = 0; i < 4; ++i)
            deriv[i] = gain_ * (1 - tanhus[i] * tanhus[i]);
        // derivation of  tanh(x - k*(g*u+s))
        T derivinput;
        if constexpr (FT == kLowPass)
            derivinput = -res_ * gain_ * (1 - tanhinput * tanhinput);
        else if (FT == kHighPass)
            derivinput = -res_ * (1 - tanhinput * tanhinput);

        // jacobian matrix coefficient
        auto a = 1 + deriv[0], b = 1 + deriv[1], c = 1 + deriv[2],
             d  = 1 + deriv[3];
        auto ap = deriv[0], bp = deriv[1], cp = deriv[2], dp = derivinput;
        if constexpr (FT == kHighPass) {
            ap = bp = cp = 1;
        }

        // compute the value of u - f(u)
        T fu[4];
        for (int i = 0; i < 4; ++i) {
            T input;
            if constexpr (FT == kLowPass)
                input = i == 0 ? tanhinput : tanhus[i - 1];
            else if constexpr (FT == kHighPass)
                input = i == 0 ? tanhinput : u[i - 1];
            fu[i] = state.u_[i] - (input - tanhus[i]);
        }

        // solve the newton value by computing direct invert matrix
        // multiplication
        // should be able to use fewer multiplications
        T solve[4];
        solve[0] = b * (c * (d * fu[0] + dp * fu[3]) + cp * dp * fu[2]) +
                   bp * cp * dp * fu[1];
        solve[1] = c * (d * (a * fu[1] + ap * fu[0]) + dp * ap * fu[3]) +
                   cp * ap * dp * fu[2];
        solve[2] = d * (a * (b * fu[2] + bp * fu[1]) + ap * bp * fu[0]) +
                   dp * ap * bp * fu[3];
        solve[3] = a * (b * (c * fu[3] + cp * fu[2]) + bp * cp * fu[1]) +
                   ap * bp * cp * fu[4];

        auto div = a * b * c * d - ap * bp * cp * dp;
        for (int i = 0; i < 4; ++i) {
            state.u_[i] -= solve[i] / div;
        }
    }

    template <class Ctxt> void process(Ctxt ctxt, State &state) const
    {
        solveLP(ctxt, state);

        T y;
        for (int i = 0; i < 4; ++i) {
            auto v      = gain_ * state.u_[i];
            y           = v + state.s_[i];
            state.s_[i] = y + v;
        }

        ctxt.setOutput(y);

        if constexpr (FT == kLowPass) ctxt.setOutput(y);
        else
            ctxt.setOutput(state.u_[3]);
    }

    PROCESSBLOCK_

  private:
    T gain_;
    T res_;
};
} // namespace dsp::va
