#pragma once

#include "Context.h"
#include "FastMath.h"
#include "Signal.h"
#include <cmath>

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

template <size_t N> static constexpr fData<N> warpGain(fData<N> freq)
{
    arrayFor(freq, i)
    {
        freq[i] = tanf(dsp::constants<float>::pi * freq[i] * 0.5f);
    }
    return freq;
}

template <size_t N, FilterType FT = kLowPass> class OnePole
{
  public:
    struct State : public fData<N> {
        State() : fData<N>{} {}
    };

    OnePole(fData<N> freq)
    {
        static_assert(FT == kLowPass || FT == kHighPass || FT == kAllPass);
        setFreq(freq);
    }
    OnePole() = default;

    void setFreq(fData<N> freq)
    {
        auto gain = warpGain(freq);
        arrayFor(freq, i) { gain_[i] = gain[i] / (1.f + gain[i]); }
    }

    template <class Ctxt> void process(Ctxt c, State &state) const
    {
        auto &x = c.getSignal();
        auto s  = state;
        inFor(x, k, i)
        {
            auto v  = (x[k][i] - s[i]) * gain_[i];
            auto lp = v + s[i];
            s[i]    = lp + v;
            if constexpr (FT == kLowPass) {
                x[k][i] = lp;
            } else {
                auto hp = x[k][i] - lp;
                if constexpr (FT == kHighPass) {
                    x[k][i] = hp;
                } else if constexpr (FT == kAllPass) {
                    auto ap = lp - hp;
                    x[k][i] = ap;
                }
            }
        }
        state = s;
    }

    PROCESSBLOCK_

    [[nodiscard]] fData<N> getGain() const
    {
        auto gain = gain_;
        if constexpr (FT == kLowPass) {
            return gain;
        } else if constexpr (FT == kHighPass) {
            arrayFor(gain, i) gain[i] = 1.f - gain[i];
        } else if constexpr (FT == kAllPass) {
            arrayFor(gain, i) gain[i] += gain[i];
        }
        return gain;
    }

  private:
    fData<N> gain_;
};

template <size_t N, FilterType FT = kLowPass> class SVF
{
  public:
    static constexpr bool kNormaLizedBandPass =
        FT == kBandPass || FT == kAllPass || FT == kNotch;

    struct State : public std::array<fData<N>, 2> {
        State() : std::array<fData<N>, 2>{} {}
    };

    SVF(fData<N> freq, fData<N> res) { setFreq(freq, res); }
    SVF() = default;

    void setFreq(fData<N> freq, fData<N> res)
    {
        auto gain = warpGain(freq);
        arrayFor(freq, i) { gain_[i] = gain[i]; }
        setRes(res);
    }
    void setFreq(fData<N> freq)
    {
        static constexpr auto kDefaultRes = 0.70710677f; // 1/sqrt(2)
        fData<N> res;
        arrayFor(res, i) { res[i] = kDefaultRes; }
        setFreq(freq, res);
    }

    void setRes(fData<N> res)
    {
        arrayFor(res, i)
        {
            assert(res[i] >= 0.f);
            denominator_[i] =
                1.f / (1.f + gain_[i] * (2.f * res[i] + gain_[i]));
            if constexpr (FT == kHighPass) {
                gains1_[i] = 2 * res[i] + gain_[i];
            }
            if constexpr (kNormaLizedBandPass) {
                inputGain_[i] = 2.f * res[i];
                if constexpr (FT == kAllPass) {
                    inputGain_[i] *= 2.f;
                }
            }
        }
    }

    /* -3db bandwidth in octave */
    template <bool PreWarp = false> void setBandWidth(fData<N> bandwith)
    {
        static_assert(FT == kBandPass);
        fData<N> res;
        arrayFor(bandwith, i)
        {
            auto bw = bandwith[i];
            if constexpr (PreWarp) {
                auto freqbwm = gain_[i] * powf(2.f, -bw * 0.5f);
                auto freqbwp = gain_[i] * powf(2.f, bw * 0.5f);
                bw           = std::log(warpGain<1>({freqbwp})[0] /
                                        warpGain<1>({freqbwm})[0]) /
                     std::log(2.f);
            }
            res[i] = std::pow(2.f, bw * 0.5f) - std::pow(2.f, -bw * 0.5f);
            res[i] *= 0.5f;
        }
        setRes(res);
    }

    template <class Ctxt> void process(Ctxt c, State &state) const
    {
        auto &in = c.getSignal();
        auto x   = in;

        if constexpr (kNormaLizedBandPass) {
            x.fromSIMD(x.toSIMD() * inputGain_.toSIMD());
        }

        auto gain        = gain_.toSIMD();
        auto denominator = denominator_.toSIMD();
        auto s1          = state[0].toSIMD();
        auto s2          = state[1].toSIMD();
        arrayFor(x, k)
        {
            auto xk = x[k].toSIMD();
            if constexpr (FT == kHighPass) {
                auto hp = (xk - gains1_.toSIMD() * s1 - s2) * denominator;
                auto v1 = gain * hp;
                auto bp = v1 + s1;
                s1      = bp + v1; // first integrator
                auto v2 = gain * bp;
                auto lp = v2 + s2;
                s2      = lp + v2; // second integrator
                in[k].fromSIMD(hp);
            } else {
                auto bp = (gain * (xk - s2) + s1) * denominator;
                if constexpr (kNormaLizedBandPass) {
                    auto bp2 = bp + bp; // first integrator;
                    s1       = bp2 - s1;
                    auto v22 = gain * bp2; // second integrator
                    s2       = v22 + s2;
                    if constexpr (FT == kBandPass) {
                        in[k].fromSIMD(bp);
                    } else if constexpr (FT == kAllPass || FT == kNotch) {
                        in[k].fromSIMD(in[k].toSIMD() - bp);
                    }
                } else if (FT == kLowPass) {
                    auto v1 = bp - s1; // first integrator
                    s1      = bp + v1;
                    auto v2 = gain * bp; // secondintegrator
                    auto lp = v2 + s2;
                    s2      = lp + v2;
                    in[k].fromSIMD(lp);
                }
            }
        }
        state[0].fromSIMD(s1);
        state[1].fromSIMD(s2);
    }

    PROCESSBLOCK_

  private:
    fData<N> gain_;
    fData<N> denominator_;

    class Empty
    {
    };
    std::conditional_t<FT == kHighPass, fData<N>, Empty> gains1_;
    std::conditional_t<kNormaLizedBandPass, fData<N>, Empty> inputGain_;
};

template <size_t N, FilterType FT = kLowPass> class Ladder
{
  public:
    using OP    = OnePole<N, FT>;
    using State = std::array<typename OP::State, 4>;

    Ladder() = default;
    Ladder(fData<N> freq, fData<N> res) : onepole_(freq) { setFreq(freq, res); }

    void setFreq(fData<N> freq, fData<N> res = {0.f})
    {
        onepole_.setFreq(freq);
        setRes(res);
    }
    void setRes(fData<N> res)
    {
        auto gain = onepole_.getGain();
        arrayFor(res, i)
        {
            auto gainSq     = gain[i] * gain[i];
            res_[i]         = res[i];
            denominator_[i] = 1.f / (1.f + gainSq * gainSq * res[i]);
        }
    }

    template <class Ctxt> void process(Ctxt c, State &state) const
    {
        auto &x = c.getSignal();

        auto gain = onepole_.getGain();
        arrayFor(x, k)
        {
            typename Ctxt::BaseType u;
            arrayFor(x[k], i)
            {
                auto s1 = state[0][i];
                auto s2 = state[1][i];
                auto s3 = state[2][i];
                auto s4 = state[3][i];
                auto g  = gain[i];
                auto S  = s4 + g * (s3 + g * (s2 + g * s1));
                // https://www.kvraudio.com/forum/viewtopic.php?f=33&t=571909
                if constexpr (FT == kLowPass) {
                    S *= (1.f - g);
                } else if constexpr (FT == kHighPass) {
                    S *= -g;
                }
                u[i] = (x[k][i] - res_[i] * S) * denominator_[i];
            }
            for (size_t j = 0; j < 4; ++j) {
                onepole_.process(Context{&u}, state[j]);
            }
            arrayFor(x[k], i) { x[k][i] = u[i]; }
        }
    }

    PROCESSBLOCK_

  private:
    fData<N> res_;
    fData<N> denominator_;
    OP onepole_;
};
} // namespace dsp::va
