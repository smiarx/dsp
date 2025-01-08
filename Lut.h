#pragma once

#include "Signal.h"

namespace dsp
{
template <int N, int Size> class Lut
{
    /* Lookup table */
    dsp::Signal<N> table_[Size][2];

    enum
    {
        valueId = 0,
        deltaId = 1,
    };

    public:
    template<typename F>
    void fill(F func)
    {
        /* func take one float input and return Signal<N> */
        auto x0 = table_[0][valueId] = func(0);
        for(int n = 1; n < Size; ++n)
        {
            auto x1 = table_[n][valueId] = func(n/Size);
            table_[n-1][deltaId] = x1-x0;
            x0 = x1;
        }
        auto x1 = func(Size/Size);
        table_[Size-1][deltaId] = x1-x0;
    }

    template <typename Float>
    dsp::Signal<N> read(Float pos, int i) const
    {
        pos *= Size;
        auto ipos = static_cast<int> (pos);
        auto fpos = pos - ipos;

        auto x = table_[ipos][valueId][i];
        const auto& dx = table_[ipos][deltaId][i];

            x += fpos*dx;
        return x;
    };

    template <typename Float>
    dsp::Signal<N> read(Float pos) const
    {
        pos *= Size;
        auto ipos = static_cast<int> (pos);
        auto fpos = pos - ipos;

        auto x = table_[ipos][valueId];
        const auto& dx = table_[ipos][deltaId];

        for(int i = 0; i < N; ++i)
        {
            x[i] += fpos*dx[i];
        }
        return x;
    };
};
} // namespace dsp
