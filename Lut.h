#pragma once

#include "Signal.h"

namespace dsp
{
template <typename T, int Size> class Lut
{
    enum {
        valueId = 0,
        deltaId = 1,
    };

  public:
    template <typename F> void fill(F func)
    {
        /* func take one float input and return Signal<N> */
        auto x0 = table_[0][valueId] = func(0);
        for (int n = 1; n < Size; ++n) {
            auto x1 = table_[n][valueId] = func(static_cast<float>(n) / Size);
            table_[n - 1][deltaId]       = x1 - x0;
            x0                           = x1;
        }
        auto x1                   = func(1.f);
        table_[Size - 1][deltaId] = x1 - x0;
    }

    template <typename Float> auto read(Float pos, int i) const
    {
        pos *= Size;
        auto ipos = static_cast<int>(pos);
        auto fpos = pos - ipos;

        auto x         = table_[ipos][valueId][i];
        const auto &dx = table_[ipos][deltaId][i];

        x += fpos * dx;
        return x;
    };

    template <typename Float> auto read(Float pos) const
    {
        pos *= Size;
        auto ipos = static_cast<int>(pos);
        auto fpos = pos - ipos;

        auto x         = table_[ipos][valueId];
        const auto &dx = table_[ipos][deltaId];

        x += dx * fpos;
        return x;
    };

  private:
    /* Lookup table */
    T table_[Size][2] = {{0.f}};
};
} // namespace dsp
