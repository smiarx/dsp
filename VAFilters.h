#pragma once

#include "Context.h"
#include "Signal.h"
#include <cmath>

namespace dsp::va
{

/* Various implementation from
 * https://archive.org/details/the-art-of-va-filter-design-rev.-2.1.2
 */

enum FilterType {
    LowPass,
    HighPass,
    AllPass,
    BandPass,
    Notch,
};

template <int N> static constexpr Signal<N> warpGain(Signal<N> freq)
{
    arrayFor(freq, i) { freq[i] = tanf(M_PIf * freq[i] * 0.5f); }
    return freq;
}

template <int N, FilterType FT = LowPass> class OnePole
{
  public:
    using State = Signal<N>;

    OnePole(Signal<N> freq)
    {
        static_assert(FT == LowPass || FT == HighPass || FT == AllPass);
        setFreq(freq);
    }
    OnePole() = default;

    void setFreq(Signal<N> freq)
    {
        auto gain = warpGain(freq);
        arrayFor(freq, i) { gain_[i] = gain[i] / (1.f + gain[i]); }
    }

    template <class Ctxt> void process(Ctxt c, State &state) const
    {
        auto &x = c.getIn();
        auto s  = state;
        inFor(x, k, i)
        {
            auto v  = (x[k][i] - s[i]) * gain_[i];
            auto lp = v + s[i];
            s[i]    = lp + v;
            if constexpr (FT == LowPass) {
                x[k][i] = lp;
            } else {
                auto hp = x[k][i] - lp;
                if constexpr (FT == HighPass) {
                    x[k][i] = hp;
                } else if constexpr (FT == AllPass) {
                    auto ap = lp - hp;
                    x[k][i] = ap;
                }
            }
        }
        state = s;
    }

    Signal<N> getGain() const
    {
        auto gain = gain_;
        if constexpr (FT == LowPass) {
            return gain;
        } else if constexpr (FT == HighPass) {
            arrayFor(gain, i) gain[i] = 1.f - gain[i];
        } else if constexpr (FT == AllPass) {
            arrayFor(gain, i) gain[i] += gain[i];
        }
        return gain;
    }

  private:
    Signal<N> gain_;
};

template <int N, FilterType FT = LowPass> class SVF
{
  public:
    static constexpr bool NormaLizedBandPass =
        FT == BandPass || FT == AllPass || FT == Notch;
    using State = std::array<Signal<N>, 2>;

    SVF(Signal<N> freq, Signal<N> res) { setFreq(freq, res); }
    SVF() = default;

    void setFreq(Signal<N> freq, Signal<N> res)
    {
        auto gain = warpGain(freq);
        arrayFor(freq, i) { gain_[i] = gain[i]; }
        setRes(res);
    }
    void setFreq(Signal<N> freq)
    {
        static constexpr auto defaultRes = 0.70710677f; // 1/sqrt(2)
        Signal<N> res;
        arrayFor(res, i) { res[i] = defaultRes; }
        setFreq(freq, res);
    }

    void setRes(Signal<N> res)
    {
        arrayFor(res, i)
        {
            assert(res[i] > 0.f);
            denominator_[i] =
                1.f / (1.f + gain_[i] * (2.f * res[i] + gain_[i]));
            if constexpr (FT == HighPass) {
                gains1_[i] = 2 * res[i] + gain_[i];
            }
            if constexpr (NormaLizedBandPass) {
                inputGain_[i] = 2.f * res[i];
                if constexpr (FT == AllPass) {
                    inputGain_[i] *= 2.f;
                }
            }
        }
    }

    template <class Ctxt> void process(Ctxt c, State &state) const
    {
        auto &in = c.getIn();
        auto x   = in;
        auto s   = state;

        if constexpr (NormaLizedBandPass) {
            inFor(x, k, i) { x[k][i] *= inputGain_[i]; }
        }

        inFor(x, k, i)
        {
            auto s1 = s[0][i], s2 = s[1][i];
            if constexpr (FT == HighPass) {
                auto hp  = (x[k][i] - gains1_[i] * s1 - s2) * denominator_[i];
                auto v1  = gain_[i] * hp;
                auto bp  = v1 + s1;
                s1       = bp + v1; // first integrator
                auto v2  = gain_[i] * bp;
                auto lp  = v2 + s2;
                s2       = lp + v2; // second integrator
                in[k][i] = hp;
            } else {
                auto bp = (gain_[i] * (x[k][i] - s2) + s1) * denominator_[i];
                if constexpr (NormaLizedBandPass) {
                    auto bp2 = bp + bp; // first integrator;
                    s1       = bp2 - s1;
                    auto v22 = gain_[i] * bp2; // second integrator
                    s2       = v22 + s2;
                    if constexpr (FT == BandPass) {
                        in[k][i] = bp;
                    } else if constexpr (FT == AllPass || FT == Notch) {
                        in[k][i] -= bp;
                    }
                } else if (FT == LowPass) {
                    auto v1  = bp - s1; // first integrator
                    s1       = bp + v1;
                    auto v2  = gain_[i] * bp; // secondintegrator
                    auto lp  = v2 + s2;
                    s2       = lp + v2;
                    in[k][i] = lp;
                }
            }

            s[0][i] = s1;
            s[1][i] = s2;
        }
        state = s;
    }

  private:
    Signal<N> gain_;
    Signal<N> denominator_;

    class Empty
    {
    };
    typename std::conditional<FT == HighPass, Signal<N>, Empty>::type gains1_;
    typename std::conditional<NormaLizedBandPass, Signal<N>, Empty>::type
        inputGain_;
};

template <int N, FilterType FT = LowPass> class Ladder
{
  public:
    using OP    = OnePole<N, FT>;
    using State = std::array<typename OP::State, 4>;

    Ladder(Signal<N> freq, Signal<N> res) : onepole_(freq)
    {
        setFreq(freq, res);
    }

    void setFreq(Signal<N> freq, Signal<N> res = {0.f})
    {
        onepole_.setFreq(freq);
        setRes(res);
    }
    void setRes(Signal<N> res)
    {
        auto gain = onepole_.getGain();
        arrayFor(res, i)
        {
            auto gainSq     = gain[i] * gain[i];
            res_[i]         = res[i];
            denominator_[i] = 1.f / (1.f + gainSq * gainSq * res[i]);
        }
    }

    template <class Ctxt> void process(Ctxt c, State &state)
    {
        auto &x = c.getIn();

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
                if constexpr (FT == LowPass) {
                    S *= (1.f - g);
                } else if constexpr (FT == HighPass) {
                    S *= -g;
                }
                u[i] = (x[k][i] - res_[i] * S) * denominator_[i];
            }
            for (int j = 0; j < 4; ++j) {
                onepole_.process(Context{&u}, state[j]);
            }
            arrayFor(x[k], i) { x[k][i] = u[i]; }
        }
    }

  private:
    Signal<N> res_;
    Signal<N> denominator_;
    OP onepole_;
};
} // namespace dsp::va
