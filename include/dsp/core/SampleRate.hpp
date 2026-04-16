#pragma once

#include <dsp/core/Units.hpp>
#include <dsp/math/FastMath.hpp>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

// SampleRate
template <typename T> class SampleRate
{
  public:
    constexpr explicit SampleRate(freq::Hz<T> sr) :
        sr_(sr.get()), invNyquist_(T(2) / sr_),
        angularScale_(constants<T>::pi * invNyquist_)
    {
    }

    constexpr void set(freq::Hz<T> sr) noexcept { *this = SampleRate(sr); }
    [[nodiscard]] constexpr freq::Hz<T> get() const noexcept
    {
        return freq::Hz{sr_};
    }

    template <typename U>
    constexpr auto hzToAngular(const freq::Hz<U> &freq) const noexcept
    {
        return freq::Angular<U>{angularScale_ * dsp::load(freq.get())};
    }

    template <typename U>
    constexpr auto hzToNyquistRatio(const freq::Hz<U> &freq) const noexcept
    {
        return freq::NyquistRatio<U>{invNyquist_ * dsp::load(freq.get())};
    }

    template <typename U>
    constexpr auto angularToHz(const freq::Angular<U> &angular) const noexcept
    {
        return freq::Hz{dsp::load(angular.get()) / constants<U>::kTwoPi * sr_};
    }

    template <typename U>
    constexpr auto
    nyquistRatioToHz(const freq::NyquistRatio<U> &nR) const noexcept
    {
        return freq::Hz{dsp::load(nR.get()) / T(2) * sr_};
    }

    template <typename U>
    constexpr auto secsToSamples(const time::Seconds<U> &sec) const noexcept
    {
        return time::Samples<U>{sr_ * dsp::load(sec.get())};
    }

    template <typename U>
    constexpr auto hzToPeriod(const freq::Hz<U> &freq) const noexcept
    {
        return time::Samples<U>{sr_ / dsp::load(freq.get())};
    }

    template <typename U>
    constexpr auto periodToHz(const time::Samples<U> &period) const noexcept
    {
        return freq::Hz<U>{sr_ / dsp::load(period.get())};
    }

    constexpr auto operator*(T factor) const noexcept
    {
        return SampleRate(get() * factor);
    }
    constexpr auto operator/(T factor) const noexcept
    {
        return SampleRate(get() / factor);
    }

    constexpr auto operator/(const SampleRate &other) const noexcept
    {
        return sr_ / other.get().get();
    }

  private:
    T sr_;
    T invNyquist_;
    T angularScale_;
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
