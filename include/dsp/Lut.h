#pragma once

#include "Signal.h"
#include "Utils.h"

namespace dsp
{
template <typename T, size_t Size> class Lut
{
    enum {
        kValueId = 0,
        kDeltaId = 1,
    };

    static_assert(isPow2(Size), "Size must be power of 2");
    static constexpr auto kMask = Size - 1;

  public:
    template <typename F> void fill(F func)
    {
        /* func take one float input and return Signal<N> */
        auto x0 = table_[0][kValueId] = func(0);
        for (size_t n = 1; n < Size; ++n) {
            auto x1 = table_[n][kValueId] = func(static_cast<float>(n) / Size);
            table_[n - 1][kDeltaId]       = x1 - x0;
            x0                            = x1;
        }
        auto x1                    = func(1.f);
        table_[Size - 1][kDeltaId] = x1 - x0;
    }

    template <typename Float> auto read(Float pos, int i) const
    {
        pos *= Size;
        auto ipos = static_cast<int>(pos);
        auto fpos = pos - ipos;

        // ipos between 0 and Size
        ipos &= static_cast<int>(kMask);

        auto x         = table_[ipos][kValueId][i];
        const auto &dx = table_[ipos][kDeltaId][i];

        x += fpos * dx;
        return x;
    }

    template <typename Float> [[nodiscard]] auto read(Float pos) const
    {
        pos *= Size;
        auto ipos = static_cast<int>(pos);
        auto fpos = pos - static_cast<float>(ipos);

        // ipos between 0 and Size
        ipos &= static_cast<int>(kMask);

        auto x         = table_[ipos][kValueId];
        const auto &dx = table_[ipos][kDeltaId];

        x += dx * fpos;
        return x;
    }

  private:
    /* Lookup table */
    T table_[Size][2] = {};
};
} // namespace dsp
