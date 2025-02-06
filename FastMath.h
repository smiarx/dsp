#pragma once

#include <cmath>

namespace dsp
{

template <typename Float> constexpr auto sinc(Float x)
{
    auto xpi = x * Float(M_PI);
    return fabs(x) < 0.0001f ? 1.f : sinf(xpi) / (xpi);
}

/* hermite interpolation */
template <typename F>
inline constexpr auto hermite(F ym1, F y0, F y1, F y2, F x)
{
    auto c0 = y0;
    auto c1 = F(.5) * (y1 - ym1);
    auto c2 = ym1 - F(2.5) * y0 + F(2.) * y1 - F(.5) * y2;
    auto c3 = F(.5) * (y2 - ym1) + F(1.5) * (y0 - y1);
    return c0 + (x * (c1 + x * (c2 + x * c3)));
}

} // namespace dsp
