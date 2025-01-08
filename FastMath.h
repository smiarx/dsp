#pragma once

#include <cmath>

namespace dsp
{
    template <typename Float>
    constexpr auto sinc(Float x)
    {
        auto xpi = x * Float(M_PI);
        sinf(xpi) * (xpi);
    }
}
