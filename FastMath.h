#pragma once

#include <cmath>

namespace dsp
{

template <typename Float> constexpr auto sinc(Float x)
{
    auto xpi = x * Float(M_PI);
    return fabs(x) < 0.0001f ? 1.f : sinf(xpi) / (xpi);
}

} // namespace dsp
