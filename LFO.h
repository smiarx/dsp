#pragma once

#include "Signal.h"
#include <cmath>

namespace dsp
{

template <int N> class LFOSine
{
    /* second midified coupled-form oscillator as seen in dattorro effect design
     * part.3
     * */
  public:
    LFOSine()
    {
        for (int i = 0; i < N; ++i) {
            cos_[i] = 1.f;
        }
    }
    void setPhase(Signal<N> phase)
    {
        for (int i = 0; i < N; ++i) {
            // TODO modify to real value
            sin_[i] = sinf(2.f * M_PIf * phase[i]);
            cos_[i] = cosf(2.f * M_PIf * phase[i]);
        }
    }

    void setFreq(Signal<N> freq)
    {
        for (int i = 0; i < N; ++i) {
            a_[i] = 2.f * sinf(M_PIf / 2.f * freq[i]);
        }
    }

    Signal<N> process()
    {
        auto y = sin_;
        for (int i = 0; i < N; ++i) {
            cos_[i] -= a_[i] * sin_[i];
            sin_[i] += a_[i] * cos_[i];
        }
        return y;
    }

  private:
    Signal<N> a_   = {{0.f}};
    Signal<N> sin_ = {{0.f}};
    Signal<N> cos_;
};

/**
 * @brief Generate sin-like parabolic Low-Frequency Oscillator (LFO) using the
 * formula: y = x⋅(2 - |x|), x in [-2,2[
 *
 * Computation use Q30 fixed point representation for phase x on 32bits
 * integers. Thus phase automaticaly wraps around 2 to -2.
 */
template <int N> class LFOParabolic
{
  public:
    // Number of fractional bits for fixed-point arithmetic
    static constexpr auto Q = 30;
    // Fixed-point representation of 1 (2^30) */
    static constexpr unsigned int Unity = 1 << Q;

    /** @brief Default constructor initializes phase to zero. */
    LFOParabolic() = default;

    /** @brief Constructor with initial phase */
    LFOParabolic(iSignal<N> phase) : phase_(phase) {}

    /** @brief Sets the initial phase for the oscillator. */
    /** @param phase Initial phase to set. */
    void setPhase(iSignal<N> phase) { phase_ = phase; }

    /** @brief Sets the frequency for the oscillator. */
    /** @param freq Frequency from [0,1] where 1 is the Nyquist frequency
     */
    void setFreq(Signal<N> freq)
    {
#pragma omp simd
        for (int i = 0; i < N; ++i) {
            assert(freq_[i] >= 0.f && freq_[i] <= 1.f);
            // Scale frequency for fixed-point accumulation
            freq_[i] = freq[i] * 4 * Unity;
        }
    }

    /**
     * @brief Generates the parabolic waveform for the current phase and
     * frequency.
     * @return Signal<N> Output lfo signal.
     */
    Signal<N> process()
    {
        Signal<N> y;

#pragma omp simd
        for (int i = 0; i < N; ++i) {

            // Compute parabolic function
            int x        = phase_[i];
            int xHalf    = x >> (Q / 2);           // Shift to Q15
            int xAbsHalf = std::abs(x) >> (Q / 2); // Shift to Q15
            int iy       = 2 * x - xHalf * xAbsHalf;

            // Increment phase
            phase_[i] += freq_[i];

            // Convert fixed-point result to float
            y[i] = static_cast<float>(iy) / Unity;
        }

        return y;
    }

  private:
    /* Current phase state for each channel */
    iSignal<N> phase_{0};
    /* Frequency of each channel */
    iSignal<N> freq_{0};
};

} // namespace dsp
