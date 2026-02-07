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
    struct State {
        T s_{};
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
        auto s = load(state.s_);

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
        state.s_ = s;
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

        auto s1 = load(state[0].s_);
        auto s2 = load(state[1].s_);
        auto s3 = load(state[2].s_);
        auto s4 = load(state[3].s_);
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

template <typename T, FilterType FT = kLowPass> class SallenKey
{
  public:
    using OP    = OnePole<T, FT == kBandPass ? kLowPass : FT>;
    using State = std::array<typename OP::State, 2>;

    SallenKey() = default;
    SallenKey(const T &freq, const T &res) : onepole_(freq)
    {
        setFreq(freq, res);
    }

    void setFreq(const T &freq, const T &res = 0.f)
    {
        onepole_.setFreq(freq);
        setRes(res);
    }
    void setRes(const T &res)
    {
        auto gain    = onepole_.getGain();
        res_         = res;
        denominator_ = T(1) / (T(1) + res_ * gain * (gain - 1));
    }

    template <class Ctxt> void process(Ctxt c, State &state) const
    {
        auto x = c.getInput();

        auto G = onepole_.getGain();

        auto s1 = load(state[0].s_);
        auto s2 = load(state[1].s_);
        // auto S  = (g-1)*s1 + s2;
        //// https://www.kvraudio.com/forum/viewtopic.php?f=33&t=571909
        // if constexpr (FT == kLowPass || FT == kBandPass) {
        //     S *= (T(1) - g);
        // } else if constexpr (FT == kHighPass) {
        //     S *= -g;
        // }

        // x = (x - res_*S) * denominator_;

        // c.setOutput(x);
        // if constexpr(FT == kBandPass)
        //{
        //     onepole_.process(c, state[0]);
        //     auto v = c.getInput();
        //     onepole_.process(c, state[1]);
        //     c.setOutput(v - c.getInput());
        // }
        // else
        //{
        //     for (size_t j = 0; j < 2; ++j) {
        //         onepole_.process(c, state[j]);
        //     }
        // }

        /*
         *
         *
         *
         *               u       v       w
         *    xlp---->|+|-->|LP|---->|+|--->|LP|----------->y
         *             ^-             ^               |
         *             |              |               |
         *             -------------------|+|--|k⋅|----
         *                                 ^
         *                                 |
         *                                xbp
         */

        auto in2 = s1 - (G - T(1)) * res_ * s2;

        if constexpr (FT == kBandPass) in2 += x;

        if constexpr (FT == kLowPass || FT == kBandPass) {
            in2 *= (T(1) - G);
            //} else if constexpr (FT == kHighPass) {
            //    in2 *= -G;
        }

        if constexpr (FT == kLowPass || FT == kHighPass) in2 += G * x;
        in2 *= denominator_;

        c.setOutput(in2);
        onepole_.process(c, state[1]);
        auto out2 = c.getInput();

        auto out1 = in2 - res_ * out2;
        if constexpr (FT == kBandPass) out1 -= x;

        if constexpr (FT == kLowPass || FT == kBandPass)
            state[0].s_ = 2 * out1 - s1;
        // else if(FT==kHighPass)
        //     state[0].s_ = 2*G*out1+s1;
    }

  private:
    T res_;
    T denominator_;
    OP onepole_;
};

///////////////////////////////////////////////////////////////////////
///////////////////////////////
//                           //
//        NON LINEAR         //
//                           //
///////////////////////////////

template <typename T> struct Jacobian2x2 {
    /*
     * |a d|
     * |c b|
     */
    T a_, b_, c_, d_;

    inline auto newtonStep(std::array<T, 2> &x, const std::array<T, 2> &fx)
    {
        /*
         * 1 step newton method for a 4 stage filter given function values and
         * jacobian
         */

        // jacobian of x - f(x)
        auto a = T(1) - a_;
        auto b = T(1) - b_;
        // should be c = -j.c but computation is easier like that
        auto c = c_;
        auto d = d_;

        // v = x - f(x)
        T v[2];
        for (int i = 0; i < 2; ++i) v[i] = x[i] - fx[i];

        // solve the newton value by computing direct invert matrix
        // multiplication
        // should be able to use fewer multiplications

        T delta[2];
        // determinant
        auto det = a * b - d * c;
        delta[0] = b * v[0] + d * v[1];
        delta[1] = c * v[0] + a * v[1];

        assert(all(abs(det) > 1e-5));

        for (int i = 0; i < 2; ++i) {
            x[i] -= delta[i] / det;
        }
    }
};

template <typename T> struct SparseJacobian4x4 {
    /*
     * Jacobian structure (for F(x)):
     *
     * | dF0/dx0     0           0        dF0/dx3 |
     * | dF1/dx0  dF1/dx1        0           0    |
     * |    0     dF2/dx1     dF2/dx2        0    |
     * |    0        0        dF3/dx2     dF3/dx3 |
     *
     * Stored as:
     *  diag[i]   = dFi/dxi
     *  lower[i]  = dFi/dx(i+1)
     *  wrap      = dF0/dx3
     */

  private:
    T diag_[4];  // main diagonal
    T lower_[3]; // sub-diagonal (i,i-1)
    T wrap_;     // (0,3) element

  public:
    constexpr SparseJacobian4x4() = default;

    constexpr SparseJacobian4x4(T d0, T d1, T d2, T d3, T l1, T l2, T l3, T w) :
        diag_{d0, d1, d2, d3}, lower_{l1, l2, l3}, wrap_(w)
    {
    }

    always_inline auto newtonStep(std::array<T, 4> &x,
                                  const std::array<T, 4> &fx)
    {
        // Solve J * delta = x - F
        // then x -= delta

        // jacobian of x - f(x)
        const T a = T(1) - diag_[0];
        const T b = T(1) - diag_[1];
        const T c = T(1) - diag_[2];
        const T d = T(1) - diag_[3];

        const T ap = lower_[0];
        const T bp = lower_[1];
        const T cp = lower_[2];
        const T dp = wrap_;

        // v = x - f(x)
        std::array<T, 4> v;
        for (int i = 0; i < 4; ++i) v[i] = x[i] - fx[i];

        // solve the newton value by computing direct invert matrix
        // multiplication
        // should be able to use fewer multiplications

        // determinant
        auto det = a * b * c * d - ap * bp * cp * dp;
        std::array<T, 4> delta;
        delta[0] =
            b * c * (d * v[0] + dp * v[3]) + cp * dp * (b * v[2] + bp * v[1]);
        delta[1] =
            c * d * (a * v[1] + ap * v[0]) + dp * ap * (c * v[3] + cp * v[2]);
        delta[2] =
            d * a * (b * v[2] + bp * v[1]) + ap * bp * (d * v[0] + dp * v[3]);
        delta[3] =
            a * b * (c * v[3] + cp * v[2]) + bp * cp * (a * v[1] + ap * v[0]);

        assert(all(abs(det) > 1e-5));
        const T invDet = T(1) / det;
        for (int i = 0; i < 4; ++i) x[i] -= delta[i] * invDet;
    }
};

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
        std::array<T, 4> s_{};
        std::array<T, 4> y_{};
        std::array<T, 4> tanhy_{};
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
        auto &hp     = state.y_;
        auto &tanhhp = state.tanhy_;
        auto &s      = state.s_;

        std::array<T, 4> derivlp;
        std::array<T, 4> lp;
        for (int i = 0; i < 4; ++i) {
            auto gtanhhp = gain_ * tanhhp[i];
            // derivation lp = of g*tanh(hp) + s
            derivlp[i] = gain_ - gtanhhp * tanhhp[i];
            lp[i]      = gtanhhp + s[i];
        }

        std::array<T, 4> fhp;
        // jacobian matrix coefficient
        SparseJacobian4x4<T> jacobian;
        if constexpr (FT == kLowPass) {
            auto input = x - res_ * lp[3];
            fhp[0]     = input - lp[0];
            fhp[1]     = lp[0] - lp[1];
            fhp[2]     = lp[1] - lp[2];
            fhp[3]     = lp[2] - lp[3];

            jacobian = {
                T(0) - derivlp[0], T(0) - derivlp[1],
                T(0) - derivlp[2], T(0) - derivlp[3],
                derivlp[0],        derivlp[1],
                derivlp[2],        T(0) - res_ * derivlp[3],
            };

        } else if (FT == kHighPass) {
            auto input = x - res_ * hp[3];
            fhp[0]     = input - lp[0];
            fhp[1]     = hp[0] - lp[1];
            fhp[2]     = hp[1] - lp[2];
            fhp[3]     = hp[2] - lp[3];

            jacobian = {
                T(0) - derivlp[0],
                T(0) - derivlp[1],
                T(0) - derivlp[2],
                T(0) - derivlp[3],
                1,
                1,
                1,
                T(0) - res_,
            };

        } else if (FT == kBandPass) {
            auto input = x + res_ * lp[3];
            fhp[0]     = input - lp[0];
            fhp[1]     = hp[0] - lp[1];
            fhp[2]     = lp[1] - lp[2];
            fhp[3]     = hp[2] - lp[3];

            jacobian = {
                T(0) - derivlp[0],
                T(0) - derivlp[1],
                T(0) - derivlp[2],
                T(0) - derivlp[3],
                1,
                derivlp[1],
                1,
                res_ * derivlp[3],
            };
        }
        jacobian.newtonStep(hp, fhp);

        for (int i = 0; i < 4; ++i) {
            tanhhp[i] = dsp::fasttanh(hp[i]);
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

        solve(ctxt, state);
        solve(ctxt, state);

        std::array<T, 4> lp;
        for (int i = 0; i < 4; ++i) {
            auto v      = gain_ * state.tanhy_[i];
            lp[i]       = v + state.s_[i];
            state.s_[i] = lp[i] + v;
        }

        if constexpr (FT == kLowPass) ctxt.setOutput(lp[3]);
        else if (FT == kHighPass || FT == kBandPass)
            ctxt.setOutput(state.y_[3]);
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
        std::array<T, 4> s_{};
        // newton resolve unknown
        std::array<T, 4> y_{};

        struct Empty {
        };
        std::conditional_t<FT == kLowPass, std::array<T, 4>, Empty> tanhy_{};
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
        auto x  = ctxt.getInput();
        auto &s = state.s_;

        // jacobian matrix coefficient
        SparseJacobian4x4<T> jacobian;

        if constexpr (FT == kLowPass) {
            // unknown y is the lp outputs
            auto &lp     = state.y_;
            auto &tanhlp = state.tanhy_;

            // lowpass output & loopback value
            std::array<T, 4> loop;
            for (int i = 0; i < 4; ++i) {
                loop[i] = gain_ * tanhlp[i];
            }

            // derivation d/dx of loopbackvalue (tanh' = 1-tanh²)
            std::array<T, 4> derivloop;
            for (int i = 0; i < 4; ++i) {
                auto derivtanhlp = T(1) - tanhlp[i] * tanhlp[i];
                derivloop[i]     = gain_ * derivtanhlp;
            }

            // updated lp
            std::array<T, 4> flp;
            if constexpr (FT == kLowPass) {
                auto tanhinput = dsp::fasttanh(x - res_ * lp[3]);
                auto derivinput =
                    T(0) - gain_ * res_ * (T(1) - tanhinput * tanhinput);
                flp[0] = gain_ * tanhinput - loop[0] + s[0];
                flp[1] = loop[0] - loop[1] + s[1];
                flp[2] = loop[1] - loop[2] + s[2];
                flp[3] = loop[2] - loop[3] + s[3];

                jacobian = {
                    T(0) - derivloop[0], T(0) - derivloop[1],
                    T(0) - derivloop[2], T(0) - derivloop[3],
                    derivloop[0],        derivloop[1],
                    derivloop[2],        derivinput,
                };
            }

            jacobian.newtonStep(lp, flp);

            for (int i = 0; i < 4; ++i) tanhlp[i] = dsp::fasttanh(lp[i]);
        }

        else if (FT == kHighPass || FT == kBandPass) {
            auto &hp = state.y_;

            // lowpass output & loopback value
            std::array<T, 4> lp;
            std::array<T, 4> loop;
            for (int i = 0; i < 4; ++i) {
                lp[i]   = gain_ * hp[i] + s[i];
                loop[i] = dsp::fasttanh(lp[i]);
            }

            // derivation d/dhp of loopbackvalue (tanh' = 1-tanh²)
            std::array<T, 4> derivloop;
            for (int i = 0; i < 4; ++i)
                derivloop[i] = gain_ * (T(1) - loop[i] * loop[i]);

            std::array<T, 4> tanhhp;
            for (int i = 0; i < 4; ++i) {
                tanhhp[i] = dsp::fasttanh(hp[i]);
            }

            // updated hp
            std::array<T, 4> fhp;
            // jacobian matrix coefficient
            if (FT == kHighPass) {
                auto tanhinput  = dsp::fasttanh(x - res_ * hp[3]);
                auto derivinput = T(0) - res_ * (T(1) - tanhinput * tanhinput);
                fhp[0]          = tanhinput - loop[0];
                fhp[1]          = tanhhp[0] - loop[1];
                fhp[2]          = tanhhp[1] - loop[2];
                fhp[3]          = tanhhp[2] - loop[3];

                jacobian = {
                    T(0) - derivloop[0],          T(0) - derivloop[1],
                    T(0) - derivloop[2],          T(0) - derivloop[3],
                    T(1) - tanhhp[0] * tanhhp[0], T(1) - tanhhp[1] * tanhhp[1],
                    T(1) - tanhhp[2] * tanhhp[2], derivinput,
                };
            }

            else if (FT == kBandPass) {
                auto tanhinput  = dsp::fasttanh(x + res_ * hp[3]);
                auto derivinput = res_ * (T(1) - tanhinput * tanhinput);
                fhp[0]          = tanhinput - loop[0];
                fhp[1]          = loop[0] - loop[1];
                fhp[2]          = tanhhp[1] - loop[2];
                fhp[3]          = loop[2] - loop[3];

                jacobian = {
                    T(0) - derivloop[0], T(0) - derivloop[1],
                    T(0) - derivloop[2], T(0) - derivloop[3],
                    derivloop[0],        T(1) - tanhhp[1] * tanhhp[1],
                    derivloop[2],        derivinput,
                };
            }

            jacobian.newtonStep(hp, fhp);
        }
    }

    template <class Ctxt>
    void process(Ctxt ctxt, State &state, int nsteps = 2) const
    {
        for (int i = 0; i < nsteps; ++i) solve(ctxt, state);

        for (int i = 0; i < 4; ++i) {
            if constexpr (FT == kLowPass) {
                state.s_[i] = T(2) * state.y_[i] - state.s_[i];
            } else if (FT == kHighPass || FT == kBandPass) {
                auto v      = gain_ * state.y_[i];
                auto lp     = v + state.s_[i];
                state.s_[i] = lp + v;
            }
        }

        ctxt.setOutput(state.y_[3]);
    }

    PROCESSBLOCK_

  private:
    T gain_{};
    T res_{};
};

template <typename T, FilterType FT = kLowPass> class SallenKeyTransistor
{
    /*
     * implements 4 stage ladder using newton method
     *                                  ─────────────────
     *                      hp1         v               │
     *  x ───>┤tanh├───>┤+├───>┤g⋅+s1├─|+|──┬───> lp1   │
     *                   ^-                 │           │
     *                   │                  │           │
     *                   ├<────┤tanh├<──────┘           │
     *                   │                              │
     *                   │               │──────────────┘
     *                   v  hp2          │
     *                  ┤+├───>┤g⋅+s2├───┬───> lp2
     *                   ^-              │
     *                   │               v
     *                   └<────┤tanh├<───┘
     *
     *   hp1 = tanh(x) - tanh(g⋅hp1+s1 + res⋅(g⋅hp2+s2))
     *   hp2 = tanh(g⋅hp1+s1 + res⋅(g⋅hp2+s2)) - tanh(g⋅hp2+s2)
     */
  public:
    using State = struct State {
        std::array<T, 2> s_{};
        std::array<T, 2> y_{};
    };

    SallenKeyTransistor() = default;
    SallenKeyTransistor(const T &freq, const T &res) { setFreq(freq, res); }

    void setFreq(const T &freq, const T &res = 0.f)
    {
        gain_ = warpGain(freq);
        setRes(res);
    }
    void setRes(const T &res) { res_ = res; }

    template <class Ctxt> void solve(Ctxt ctxt, State &state) const
    {
        auto x   = ctxt.getInput();
        auto &hp = state.y_;
        auto &s  = state.s_;

        auto input = dsp::fasttanh(x);
        auto gain0 = gain_ * hp[0] + s[0];
        auto gain1 = gain_ * hp[1] + s[1];
        auto loop0 = dsp::fasttanh(gain0 + res_ * gain1);
        // auto tanhloop = dsp::fasttanh(res_*gain1);
        // auto loop0 = dsp::fasttanh(gain0+tanhloop);
        auto loop1 = dsp::fasttanh(gain1);

        auto derivloop00 = gain_ * (1 - loop0 * loop0);
        auto derivloop01 = res_ * derivloop00;
        // auto derivloop01 = res_*(1-tanhloop*tanhloop)*derivloop00;
        auto derivloop1 = gain_ * (1 - loop1 * loop1);

        std::array<T, 2> fhp = {input - loop0, loop0 - loop1};
        Jacobian2x2<T> jacobian;
        jacobian = {-derivloop00, derivloop01 - derivloop1, derivloop00,
                    -derivloop01};

        jacobian.newtonStep(hp, fhp);

        //// lowpass output & loopback value
        // T lp[2];
        // T loop[2];
        // for (int i = 0; i < 2; ++i) {
        //     lp[i]   = gain_ * hp[i] + s[i];
        //     loop[i] = dsp::fasttanh(lp[i]);
        // }

        //// derivation d/dhp of loopbackvalue (tanh' = 1-tanh²)
        // T derivloop[2];
        // for (int i = 0; i < 2; ++i)
        //     derivloop[i] = gain_ * (1 - loop[i] * loop[i]);

        //// jacobian matrix coefficient
        // Jacobian2x2<T> jacobian;

        // T fhp[2];
        // if constexpr (FT == kLowPass)
        //{
        //     auto tanhinput = dsp::fasttanh(x - res_ * hp[1]);
        //     auto derivinput = -res_ * (1 - tanhinput * tanhinput);

        //    auto input2 = dsp::fasttanh(lp[0]+res_*hp[1]);
        //    auto derivinput20 = gain_*(1-input2*input2);
        //    auto derivinput21 = res_*(1-input2*input2);

        //    // compute the updated hp
        //    fhp[0] = tanhinput - loop[0];
        //    fhp[1] = input2 - loop[1];

        //    jacobian = {-derivloop[0], derivinput21-derivloop[1],
        //    derivinput20, derivinput};
        //}
        // solveNewton(hp, jacobian, fhp);

        //// lowpass output & loopback value
        // T lp[2];
        // T loop[2];
        // for (int i = 0; i < 2; ++i) {
        //     lp[i]   = gain_ * hp[i] + s[i];
        //     loop[i] = dsp::fasttanh(lp[i]);
        // }

        //// derivation d/dhp of loopbackvalue (tanh' = 1-tanh²)
        // T derivloop[2];
        // for (int i = 0; i < 2; ++i)
        //     derivloop[i] = gain_ * (1 - loop[i] * loop[i]);

        //// jacobian matrix coefficient
        // Jacobian2x2<T> jacobian;

        // T fhp[2];
        // if constexpr (FT == kLowPass)
        //{
        //     auto tanhinput = dsp::fasttanh(x + res_ * hp[1]);
        //     auto derivinput = res_ * (1 - tanhinput * tanhinput);

        //    // compute the updated hp
        //    fhp[0] = tanhinput - loop[0];
        //    fhp[1] = loop[0] - loop[1];

        //    jacobian = {-derivloop[0], -derivloop[1], derivloop[0],
        //    derivinput};
        //}
        // solveNewton(hp, jacobian, fhp);
    }

    template <class Ctxt> void process(Ctxt ctxt, State &state) const
    {
        solve(ctxt, state);
        solve(ctxt, state);

        T y[2];
        for (int i = 0; i < 2; ++i) {
            auto v      = gain_ * state.y_[i];
            y[i]        = v + state.s_[i];
            state.s_[i] = y[i] + v;
        }

        if constexpr (FT == kLowPass) ctxt.setOutput(y[1]);
        // else if (FT == kHighPass || FT == kBandPass)
        //     ctxt.setOutput(state.y_[3]);
    }

    PROCESSBLOCK_

  private:
    T gain_;
    T res_;
};

template <typename T, FilterType FT = kLowPass> class SallenKeyDiode
{
    /*
     * implements 4 stage ladder using newton method
     *                         ┌───────────────┐
     *             hp1         v               │
     *  x ────>┤+├───>┤g⋅+s1├─|+|──┬───> lp1   │
     *          ^-                 │           │
     *          │                  v           △⋅res
     *          └<────┤sinh├<──────┤           │
     *                             │         │tanh│
     *          ┌──────────────────┘           ^
     *          │                              │
     *          │                              │
     *          │               ┌──────────────┘
     *          v  hp2          │
     *         ┤+├───>┤g⋅+s2├──>┼───> lp2
     *          ^-              │
     *          │               v
     *          └─────┤sinh├<───┘
     *
     *   lp1 = g⋅(x - sinh(lp1))  + s + res⋅tanh(lp2))
     *   lp2 = g⋅(lp1 - sinh(lp2)) + s
     */
  public:
    using State = struct State {
        std::array<T, 2> s_{};
        std::array<T, 2> y_{};
        std::array<T, 2> sinhy_{};
    };

    static T fsinh(T x)
    {
        auto x2 = x * x;
        return x * (T(1) + x2 * (T(.16666666666666666) +
                                 x2 * (T(0.008333333333333333) +
                                       x2 * T(0.001851851851851852))));
    }
    static T dfsinh_dx(T x)
    {
        auto x2 = x * x;
        return T(1) + x2 * (T(.5) + x2 * (T(0.041666666666666664) +
                                          x2 * T(0.001388888888888889)));
    }

    SallenKeyDiode() = default;
    SallenKeyDiode(const T &freq, const T &res) { setFreq(freq, res); }

    void setFreq(const T &freq, const T &res = 0.f)
    {
        gain_ = warpGain(freq);
        setRes(res);
    }
    void setRes(const T &res) { res_ = res; }

    template <class Ctxt> void solve(Ctxt ctxt, State &state) const
    {
        auto x       = ctxt.getInput();
        auto &lp     = state.y_;
        auto &sinhlp = state.sinhy_;
        auto &s      = state.s_;

        auto loop = fasttanh(lp[1]);

        auto derivloop0 = gain_ * dfsinh_dx(lp[0]);
        auto derivloop1 = gain_ * dfsinh_dx(lp[1]);
        auto derivloop  = res_ * (T(1) - loop * loop);

        T xlp, xhp;
        if constexpr (FT == kLowPass) {
            xlp = x, xhp = 0;
        } else if constexpr (FT == kBandPass || FT == kHighPass) {
            xlp = 0, xhp = x;
        }

        std::array<T, 2> flp = {gain_ * (xlp - sinhlp[0]) + s[0] +
                                    res_ * (xhp + loop),
                                gain_ * (lp[0] - sinhlp[1]) + s[1]};
        Jacobian2x2<T> jacobian;
        jacobian = {-derivloop0, -derivloop1, gain_, derivloop};

        jacobian.newtonStep(lp, flp);

        sinhlp[0] = fsinh(lp[0]);
        sinhlp[1] = fsinh(lp[1]);
    }

    template <class Ctxt> void process(Ctxt ctxt, State &state) const
    {
        solve(ctxt, state);
        solve(ctxt, state);

        auto x = ctxt.getInput();
        if constexpr (FT == kLowPass)
            state.s_[0] += T(2) * gain_ * (x - state.sinhy_[0]);
        else if (FT == kHighPass || FT == kBandPass)
            state.s_[0] -= T(2) * gain_ * (state.sinhy_[0]);

        state.s_[1] = T(2) * state.y_[1] - state.s_[1];

        if constexpr (FT == kLowPass || FT == kBandPass)
            ctxt.setOutput(state.y_[1]);
        else if constexpr (FT == kHighPass)
            ctxt.setOutput(state.y_[0] - state.y_[1]);
    }

    PROCESSBLOCK_

  private:
    T gain_;
    T res_;
};

template <typename T, FilterType FT = kLowPass> class SallenKeyDiode2
{
    /*
     * implements 4 stage ladder using newton method
     *                         ┌───────────────┐
     *             hp1         v               │
     *  x ────>┤+├───>┤g⋅+s1├─|+|──┬───> lp1   │
     *          ^-                 │           │
     *          │                  v           △⋅res
     *          └<───┤tanh⁻¹├<─────┤           │
     *                             │           │
     *          ┌──────────────────┘           ^
     *          │                              │
     *          │                              │
     *          │               ┌──────────────┘
     *          v  hp2          │
     *         ┤+├───>┤g⋅+s2├──>┼───> lp2
     *          ^-              │
     *          │               v
     *          └───┤tanh⁻¹├<───┘
     *
     * lp1 = tanh(x-(lp1-res⋅lp2-s1)/g)
     * lp2 = tanh(lp1 - (lp2-s2)/g)
     */
  public:
    using State = struct State {
        std::array<T, 2> s_{};
        std::array<T, 2> y_{};
    };

    SallenKeyDiode2() = default;
    SallenKeyDiode2(const T &freq, const T &res) { setFreq(freq, res); }

    void setFreq(const T &freq, const T &res = 0.f)
    {
        invgain_ = T(1) / warpGain(freq);
        setRes(res);
    }
    void setRes(const T &res) { res_ = res; }

    template <class Ctxt> void solve(Ctxt ctxt, State &state) const
    {
        auto x   = ctxt.getInput();
        auto &lp = state.y_;
        auto &s  = state.s_;

        /*
         * lp1 = tanh(x-(lp1-res⋅lp2-s1)/g)
         * lp2 = tanh(lp1 - (lp2-s2)/g)
         */
        auto flp0 = dsp::tanh(x - (lp[0] - res_ * lp[1] - s[0]) * invgain_);
        auto flp1 = dsp::tanh(lp[0] - (lp[1] - s[1]) * invgain_);

        auto d0 = (T(1) - flp0 * flp0);
        auto d1 = (T(1) - flp1 * flp1);

        std::array<T, 2> flp = {flp0, flp1};
        Jacobian2x2<T> jacobian;
        jacobian = {-invgain_ * d0, -invgain_ * d1, d0, res_ * invgain_ * d1};

        jacobian.newtonStep(lp, flp);
    }

    template <class Ctxt> void process(Ctxt ctxt, State &state) const
    {
        solve(ctxt, state);
        solve(ctxt, state);
        solve(ctxt, state);
        solve(ctxt, state);
        solve(ctxt, state);
        solve(ctxt, state);

        state.s_[0] = T(2) * (state.y_[0] - res_ * state.y_[1]) - state.s_[0];
        state.s_[1] = T(2) * state.y_[1] - state.s_[1];

        ctxt.setOutput(state.y_[1]);
    }

    PROCESSBLOCK_

  private:
    T invgain_;
    T res_;
};
} // namespace dsp::va
