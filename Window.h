#pragma once

#include "FastMath.h"

namespace dsp
{
namespace window
{

template <int Attenuation> class Kaiser
{
  public:
    Kaiser() = delete;
    template <typename F> static constexpr auto generate(F x)
    {
        constexpr F beta  = computeBeta<F>(Attenuation);
        constexpr F denom = F(1) / zerothOrderBessel(beta);
        F K               = beta * sqrtf(F(1) - x * x);
        F num             = zerothOrderBessel(K);

        return num * denom;
    }

    template <typename F> static constexpr F computeBeta(int Att)
    {
        if (Att < 21) {
            return 0;
        } else if (Att < 50) {
            return pow(0.5842 * (Att - 21), 0.4) + 0.07886 * (Att - 21);
        } else {
            return 0.1102 * (Att - 8.7);
        }
    }
};

} // namespace window
} // namespace dsp
