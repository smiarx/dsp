#pragma once

#include "FastMath.h"
#include <cassert>

namespace dsp::lfo
{
using namespace loadfuncs;

template <typename T> class Sine
{
    /* second modified coupled-form oscillator as seen in dattorro effect design
     * part.3
     * */
  public:
    Sine() = default;
    Sine(const T &freq) { setFreq(freq); }
    Sine(const T &freq, const T &phase) { setFreq(freq, phase); }

    void setFreq(const T &freq)
    {
        auto w_2 = dsp::constants<T>::pi / 2 * load(freq);
        a_       = dsp::sin(w_2) * 2;
    }

    void setFreq(const T &freq, const T &phase)
    {
        setFreq(freq);

        auto w_2   = dsp::constants<T>::pi / 2 * load(freq);
        auto scale = T(1) / dsp::cos(w_2);

        auto wphase = load(phase) * constants<T>::pi * 2;
        sin_        = scale * dsp::sin(wphase);
        cos_        = scale * dsp::cos(wphase - w_2);
    }

    [[nodiscard]] auto process()
    {
        auto y = load(sin_);
        cos_ -= load(a_) * load(sin_);
        sin_ += load(a_) * load(cos_);
        return y;
    }

  private:
    T a_{};
    T sin_{};
    T cos_{1.};
};

/**
 * @brief Generate sin-like parabolic Low-Frequency Oscillator (LFO) using the
 * formula: y = x⋅(2 - |x|), x in [-2,2[
 *
 * Computation use Q30 fixed point representation for phase x on 32bits
 * integers. Thus phase automaticaly wraps around 2 to -2.
 */
template <typename T> class Parabolic
{
  public:
    using iT = intType<T>;
    // Number of fractional bits for fixed-point arithmetic
    static constexpr auto kQ = 30;
    // Fixed-point representation of 1 (2^30) */
    static constexpr unsigned int kUnity = 1 << kQ;

    /** @brief Default constructor initializes phase to zero. */
    Parabolic() = default;

    /** @brief Constructor with initial phase */
    Parabolic(const T &freq) { setFreq(freq); }
    Parabolic(const T &freq, const iT &phase) : phase_(phase) { setFreq(freq); }

    /** @brief Sets the initial phase for the oscillator. */
    /** @param phase Initial phase to set. */
    void setPhase(const iT &phase) { phase_ = phase; }

    /** @brief Sets the frequency for the oscillator. */
    /** @param freq Frequency from [0,1] where 1 is the Nyquist frequency
     */
    void setFreq(const T &freq)
    {
        auto f = load(freq);
        assert(all(f >= 0 && f <= 1));
        freq_ = toInt(f * 4 * kUnity);
    }

    /**
     * @brief Generates the parabolic waveform for the current phase and
     * frequency.
     * @return Signal<N> Output lfo signal.
     */
    auto process()
    {
        auto x        = load(phase_);
        auto xHalf    = x >> (kQ / 2);      // Shift to Q15
        auto xAbsHalf = abs(x) >> (kQ / 2); // Shift to Q15
        auto iy       = 2 * x - xHalf * xAbsHalf;

        phase_ += load(freq_);
        auto y = toFloat<baseType<T>>(iy) / kUnity;

        return y;
    }

  private:
    /* Current phase state for each channel */
    iT phase_{};
    /* Frequency of each channel */
    iT freq_{};
};

} // namespace dsp::lfo
