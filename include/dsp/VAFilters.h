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

template <typename T> struct SparseJacobian {
    T a, b, c, d, ap, bp, cp, dp;
};

template <typename T>
static inline auto solveNewton4Stages(T *x, const SparseJacobian<T> &j,
                                      const T *fx)
{
    /*
     * 1 step newton method for a 4 stage filter whose jacobian looks like
     * |a  0  0  d'|
     * |a' b  0  0 |
     * |0  b' c  0 |
     * |0  0  c' d |
     */

    // jacobian of x - f(x)
    auto a = 1 - j.a;
    auto b = 1 - j.b;
    auto c = 1 - j.c;
    auto d = 1 - j.d;
    // should be ap = -j.ap but computation is easier like that
    auto ap = j.ap;
    auto bp = j.bp;
    auto cp = j.cp;
    auto dp = j.dp;

    T solve[4];

    // v = x - f(x)
    T v[4];
    for (int i = 0; i < 4; ++i) v[i] = x[i] - fx[i];

    // solve the newton value by computing direct invert matrix
    // multiplication
    // should be able to use fewer multiplications
    auto div = a * b * c * d - ap * bp * cp * dp;
    solve[0] =
        b * (c * (d * v[0] + dp * v[3]) + cp * dp * v[2]) + bp * cp * dp * v[1];
    solve[1] =
        c * (d * (a * v[1] + ap * v[0]) + dp * ap * v[3]) + cp * ap * dp * v[2];
    solve[2] =
        d * (a * (b * v[2] + bp * v[1]) + ap * bp * v[0]) + dp * ap * bp * v[3];
    solve[3] =
        a * (b * (c * v[3] + cp * v[2]) + bp * cp * v[1]) + ap * bp * cp * v[0];

    for (int i = 0; i < 4; ++i) {
        x[i] -= solve[i] / div;
    }
}

template <typename T, FilterType FT = kLowPass> class LadderOTA
{
    /*
     * implements 4 stage ladder using newton method
     *
     *                   hp1
     *  x ───>┤+├───>│+├───>┤tanh├──>┤g⋅+s1├─┬──> lp1
     *         ^-     ^-                     │
     *         │      │                      │
     *         │      ├<─────────────────────┘
     *         │      │
     *         │      v  hp2
     *         │     │+├───>┤tanh├──>┤g⋅+s2├─┬──> lp2
     *         │      ^-                     │
     *         │      │                      │
     *         △*res  ├<─────────────────────┘
     *         │      │
     *         │      v  hp3
     *         │     │+├───>┤tanh├──>┤g⋅+s3├─┬──> lp3
     *         │      ^-                     │
     *         │      │                      │
     *         │      ├<─────────────────────┘
     *         │      │
     *         │      v  hp4
     *         │     │+├───>┤tanh├──>┤g⋅+s4├─┬──> lp4
     *         │      ^-                     │
     *         │      │                      │
     *         └──────┴<─────────────────────┘
     *
     * highpass and bandpass are computed reordering the the inputs of the
     * stages (hp->hp->hp->hp and lp->hp->lp->hp)
     *
     */
  public:
    using State = struct State {
        T s_[4]{};
        T hp_[4]{};
        T tanhhp_[4]{};
    };

    LadderOTA() = default;
    LadderOTA(const T &freq, const T &res) { setFreq(freq, res); }

    void setFreq(const T &freq, const T &res = 0.f)
    {
        gain_ = warpGain(freq);
        setRes(res);
    }
    void setRes(const T &res) { res_ = res; }

    template <class Ctxt> void solve(Ctxt ctxt, State &state) const
    {
        auto x       = ctxt.getInput();
        auto &hp     = state.hp_;
        auto &tanhhp = state.tanhhp_;
        auto &s      = state.s_;

        T derivlp[4];
        T lp[4];
        for (int i = 0; i < 4; ++i) {
            auto gtanhhp = gain_ * tanhhp[i];
            // derivation lp = of g*tanh(hp) + s
            derivlp[i] = gain_ - gtanhhp * tanhhp[i];
            lp[i]      = gtanhhp + s[i];
        }

        T fhp[4];
        for (int i = 0; i < 4; ++i) {
            T input;
            if constexpr (FT == kLowPass) {
                input = i == 0 ? x - res_ * lp[3] : lp[i - 1];
            } else if (FT == kHighPass) {
                input = i == 0 ? x - res_ * hp[3] : hp[i - 1];
            } else if (FT == kBandPass) {
                input = i == 0             ? x + res_ * lp[3]
                        : (i - 1) % 2 == 0 ? hp[i - 1]
                                           : lp[i - 1];
            }
            fhp[i] = input - lp[i];
        }

        // jacobian matrix coefficient
        SparseJacobian<T> jacobian;
        if constexpr (FT == kLowPass) {
            jacobian = {
                -derivlp[0], -derivlp[1], -derivlp[2], -derivlp[3],
                derivlp[0],  derivlp[1],  derivlp[2],  -res_ * derivlp[3],
            };
        } else if (FT == kHighPass) {
            jacobian = {
                -derivlp[0], -derivlp[1], -derivlp[2], -derivlp[3],
                1,           1,           1,           -res_,
            };
        } else if (FT == kBandPass) {
            jacobian = {
                -derivlp[0], -derivlp[1], -derivlp[2], -derivlp[3],
                1,           derivlp[1],  1,           res_ * derivlp[3],
            };
        }
        solveNewton4Stages(hp, jacobian, fhp);

        for (int i = 0; i < 4; ++i) {
            tanhhp[i] = dsp::fasttanh(hp[i]);
        }
    }

    template <class Ctxt> void process(Ctxt ctxt, State &state) const
    {
        solve(ctxt, state);

        T lp[4];
        for (int i = 0; i < 4; ++i) {
            auto v      = gain_ * state.tanhhp_[i];
            lp[i]       = v + state.s_[i];
            state.s_[i] = lp[i] + v;
        }

        if constexpr (FT == kLowPass) ctxt.setOutput(lp[3]);
        else if (FT == kHighPass || FT == kBandPass)
            ctxt.setOutput(state.hp_[3]);
    }

    PROCESSBLOCK_

  private:
    T gain_;
    T res_;
};

template <typename T, FilterType FT = kLowPass> class LadderTransistor
{
    /*
     * implements 4 stage ladder using newton method
     *
     *                             hp1
     *  x ───>┤+├───>┤tanh├───>┤+├───>┤g⋅+s1├───┬───> lp1
     *         ^-               ^-              │
     *         │                │               │
     *         │                ├<────┤tanh├<───┘
     *         │                │
     *         │                v  hp2
     *         │               │+├───>┤g⋅+s2├───┬───> lp2
     *         │                ^-              │
     *         │                │               │
     *         △*res            ├<────┤tanh├<───┘
     *         │                │
     *         │                v  hp3
     *         │               │+├───>┤g⋅+s3├───┬───> lp3
     *         │                ^-              │
     *         │                │               │
     *         │                ├<────┤tanh├<───┘
     *         │                │
     *         │                v  hp3
     *         │               │+├───>┤g⋅+s4├───┬───> lp4
     *         │                ^-              │
     *         │                │               v
     *         │                └<────┤tanh├<───┤
     *         │                                │
     *         └────────────────────────────────┘
     *
     * highpass and bandpass are computed reordering the the inputs of the
     * stages (hp->hp->hp->hp and lp->hp->lp->hp, in these case a tanh need to
     * be added for the input, becaused it has been removed in the lowpass graph
     * because of redundancy between loop and input)
     */
  public:
    using State = struct State {
        T s_[4]{};
        T hp_[4]{};
    };

    LadderTransistor() = default;
    LadderTransistor(const T &freq, const T &res) { setFreq(freq, res); }

    void setFreq(const T &freq, const T &res = 0.f)
    {
        gain_ = warpGain(freq);
        setRes(res);
    }
    void setRes(const T &res) { res_ = res; }

    template <class Ctxt> void solve(Ctxt ctxt, State &state) const
    {
        auto x   = ctxt.getInput();
        auto &hp = state.hp_;
        auto &s  = state.s_;

        // lowpass output & loopback value
        T lp[4];
        T loop[4];
        for (int i = 0; i < 4; ++i) {
            lp[i]   = gain_ * hp[i] + s[i];
            loop[i] = dsp::fasttanh(lp[i]);
        }

        // derivation d/dhp of loopbackvalue (tanh' = 1-tanh²)
        T derivloop[4];
        for (int i = 0; i < 4; ++i)
            derivloop[i] = gain_ * (1 - loop[i] * loop[i]);

        // filter input
        T tanhinput;
        if constexpr (FT == kLowPass)
            tanhinput = dsp::fasttanh(x - res_ * lp[3]);
        else if (FT == kHighPass)
            tanhinput = dsp::fasttanh(x - res_ * hp[3]);
        else if (FT == kBandPass)
            tanhinput = dsp::fasttanh(x + res_ * hp[3]);

        // derivation of filter input
        T derivinput;
        if constexpr (FT == kLowPass)
            derivinput = -res_ * gain_ * (1 - tanhinput * tanhinput);
        else if (FT == kHighPass)
            derivinput = -res_ * (1 - tanhinput * tanhinput);
        else if (FT == kBandPass)
            derivinput = res_ * (1 - tanhinput * tanhinput);

        [[maybe_unused]] T tanhhp[4];
        if (FT == kHighPass || FT == kBandPass) {
            for (int i = 0; i < 4; ++i) {
                tanhhp[i] = dsp::fasttanh(hp[i]);
            }
        }

        // compute the updated hp
        T fhp[4];
        for (int i = 0; i < 4; ++i) {
            T input;
            if constexpr (FT == kLowPass)
                input = i == 0 ? tanhinput : loop[i - 1];
            else if constexpr (FT == kHighPass)
                input = i == 0 ? tanhinput : tanhhp[i - 1];
            else if constexpr (FT == kBandPass)
                input = i == 0             ? tanhinput
                        : (i - 1) % 2 == 0 ? loop[i - 1]
                                           : tanhhp[i - 1];
            fhp[i] = input - loop[i];
        }

        // jacobian matrix coefficient
        SparseJacobian<T> jacobian;
        if (FT == kLowPass) {
            jacobian = {
                -derivloop[0], -derivloop[1], -derivloop[2], -derivloop[3],
                derivloop[0],  derivloop[1],  derivloop[2],  derivinput,
            };
        } else if constexpr (FT == kHighPass) {
            jacobian = {
                -derivloop[0],
                -derivloop[1],
                -derivloop[2],
                -derivloop[3],
                1 - tanhhp[0] * tanhhp[0],
                1 - tanhhp[1] * tanhhp[1],
                1 - tanhhp[2] * tanhhp[2],
                derivinput,
            };
        } else if constexpr (FT == kBandPass) {
            jacobian = {
                -derivloop[0], -derivloop[1], -derivloop[2],
                -derivloop[3], derivloop[0],  1 - tanhhp[1] * tanhhp[1],
                derivloop[2],  derivinput,
            };
        }
        solveNewton4Stages(hp, jacobian, fhp);
    }

    template <class Ctxt> void process(Ctxt ctxt, State &state) const
    {
        solve(ctxt, state);

        T y[4];
        for (int i = 0; i < 4; ++i) {
            auto v      = gain_ * state.hp_[i];
            y[i]        = v + state.s_[i];
            state.s_[i] = y[i] + v;
        }

        if constexpr (FT == kLowPass) ctxt.setOutput(y[3]);
        else if (FT == kHighPass || FT == kBandPass)
            ctxt.setOutput(state.hp_[3]);
    }

    PROCESSBLOCK_

  private:
    T gain_;
    T res_;
};
} // namespace dsp::va
