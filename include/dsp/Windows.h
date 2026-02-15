#pragma once

#include "FastMath.h"

namespace dsp::windows
{

template <int Attenuation> class Kaiser
{
  public:
    Kaiser() = delete;
    template <typename F> static constexpr auto generate(F x)
    {
        constexpr F kBeta  = computeBeta<F>(Attenuation);
        constexpr F kDenom = F(1) / zerothOrderBessel(kBeta);
        F K                = kBeta * dsp::sqrt(F(1) - x * x);
        F num              = zerothOrderBessel(K);

        return num * kDenom;
    }

    template <typename F> static constexpr F computeBeta(int att)
    {
        auto fAtt = static_cast<F>(att);
        if (att < 21) {
            return F(0);
        } else if (att < 50) {
            return std::pow(F(0.5842) * (fAtt - F(21)), F(0.4)) +
                   F(0.07886) * (fAtt - F(21));
        } else {
            return F(0.1102) * (fAtt - F(8.7));
        }
    }
};

class Hann
{
  public:
    Hann() = delete;
    template <typename F> static constexpr auto generate(F x)
    {
        return F(0.5) + F(0.5) * std::cos(constants<F>::pi * x);
    }
};

class Hamming
{
  public:
    Hamming() = delete;
    template <typename F> static constexpr auto generate(F x)
    {
        // w(n)=0.54−0.46cos(2πn/N−1)
        return F(0.54) - F(0.46) * cos(constants<F>::pi * (x + 1));
    }
};

} // namespace dsp::windows
