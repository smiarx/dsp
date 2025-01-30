#pragma once

#include "Signal.h"

namespace dsp
{

template <int N> class LFOParabolic
{

    /*
     * Generate sin-like lfo with two parabolas
     *      | y =  sign(x) ⋅ [2⋅|x| - x²]
     *  x in [-2,2]
     */
  public:
    static constexpr auto Q             = 30;
    static constexpr unsigned int Unity = 1 << Q;
    static constexpr auto Lim           = Unity << 1;
    static constexpr auto MaskLim       = Lim - 1;

    LFOParabolic() = default;
    LFOParabolic(iSignal<N> phase) : phase_(phase) {}

    void setFreq(Signal<N> freq)
    {
        for (int i = 0; i < N; ++i) {
            freq_[i] = freq[i] * 4 * Unity;
        }
    }

    Signal<N> process()
    {
        Signal<N> y;

#pragma omp simd
        for (int i = 0; i < N; ++i) {
            int x  = phase_[i] & MaskLim;
            int x2 = x >> (Q / 2);
            int iy = 2 * x - (x2 * x2);

            if (phase_[i] < 0) iy = -iy;

            phase_[i] += freq_[i];

            y[i] = static_cast<float>(iy) / Unity;
        }
        return y;
    }

  private:
    iSignal<N> phase_{0};
    iSignal<N> freq_{0};
};

} // namespace dsp
