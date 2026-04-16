#pragma once

#include <dsp/math/FastMath.hpp>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <template <typename> class Impl, typename T> class Unit
{
  public:
    using type = T;
    constexpr explicit Unit(T value) : value_{value} {}
    [[nodiscard]] constexpr auto get() const noexcept { return value_; }

    ////// operators
    template <typename U> constexpr auto operator*(U factor) const
    {
        auto r = value_ * factor;
        return Impl<decltype(r)>{r};
    }
    template <typename U> constexpr auto &operator*=(U factor)
    {
        value_ *= factor;
        return *this;
    }

    template <typename U> constexpr auto operator/(U factor) const
    {
        auto r = value_ / factor;
        return Impl<decltype(r)>{r};
    }

    template <typename U> constexpr auto operator/(Impl<U> other) const
    {
        return value_ / other.get();
    }

    template <typename U> constexpr auto operator+(Impl<U> other) const
    {
        auto r = value_ + other.get();
        return Impl<decltype(r)>{r};
    }
    template <typename U> constexpr auto &operator+=(Impl<U> other)
    {
        value_ += other.get();
        return *this;
    }

    template <typename U> constexpr auto operator-(Impl<U> other) const
    {
        auto r = value_ - other.get();
        return Impl<decltype(r)>{r};
    }

    template <typename U> constexpr bool operator==(Impl<U> other) const
    {
        return value_ == other.get();
    }

    template <typename U> constexpr bool operator!=(Impl<U> other) const
    {
        return value_ != other.get();
    }

    template <typename U> constexpr bool operator<(Impl<U> other) const
    {
        return value_ < other.get();
    }

    template <typename U> constexpr bool operator<=(Impl<U> other) const
    {
        return value_ <= other.get();
    }

    template <typename U> constexpr bool operator>(Impl<U> other) const
    {
        return value_ > other.get();
    }

    template <typename U> constexpr bool operator>=(Impl<U> other) const
    {
        return value_ >= other.get();
    }

  private:
    T value_;
};

template <template <typename> class Impl, typename T, typename U>
constexpr auto operator*(U factor, Unit<Impl, T> unit)
{
    auto r = factor * unit.get();
    return Impl<decltype(r)>{r};
}

template <template <typename> class Impl, typename T, typename U>
constexpr auto min(Unit<Impl, T> unit1, Unit<Impl, U> unit2)
{
    auto r = dsp::min(dsp::load(unit1.get()), dsp::load(unit2.get()));
    return Impl<decltype(r)>{r};
}

template <template <typename> class Impl, typename T, typename U>
constexpr auto max(Unit<Impl, T> unit1, Unit<Impl, U> unit2)
{
    auto r = dsp::max(dsp::load(unit1.get()), dsp::load(unit2.get()));
    return Impl<decltype(r)>{r};
}

template <template <typename> class Impl, typename T>
constexpr auto abs(Unit<Impl, T> unit)
{
    auto r = dsp::abs(dsp::load(unit.get()));
    return Impl<decltype(r)>{r};
}

/////////////////////////////////////////////////////////////////////

namespace time
{
template <typename T> struct Samples;
}

namespace freq
{

template <typename T> struct Hz : public Unit<Hz, T> {
    constexpr Hz(T value) : Unit<Hz, T>{value} {}
};

template <typename T> struct Angular;
template <typename T> struct NyquistRatio;

template <typename T> struct Angular : public Unit<Angular, T> {
    constexpr explicit Angular(T value) : Unit<Angular, T>{value} {}
    constexpr explicit Angular(NyquistRatio<T> ratio);
};

template <typename T> struct NyquistRatio : public Unit<NyquistRatio, T> {
    constexpr explicit NyquistRatio(T value) : Unit<NyquistRatio, T>{value} {}
    constexpr explicit NyquistRatio(Angular<T> ang);

    [[nodiscard]] constexpr dsp::time::Samples<T> toPeriod() const;
};

template <typename T>
constexpr Angular<T>::Angular(NyquistRatio<T> ratio) :
    Angular(constants<T>::pi * dsp::load(ratio.get()))
{
}

template <typename T>
constexpr NyquistRatio<T>::NyquistRatio(Angular<T> ang) :
    NyquistRatio(constants<T>::invPi * dsp::load(ang.get()))
{
}

template <typename T> [[nodiscard]] auto warp(Angular<T> angular) noexcept
{
    return dsp::tan(T(0.5) * load(angular.get()));
}

template <typename T> [[nodiscard]] auto warp(NyquistRatio<T> ratio) noexcept
{
    return dsp::tan(constants<T>::pi_2 * load(ratio.get()));
}

template <class F>
concept SampleRateIndependent = requires(F freq) { dsp::freq::warp(freq); };

template <class F>
concept Frequency = std::is_same_v<Hz<typename F::type>, F> or
                    std::is_same_v<NyquistRatio<typename F::type>, F> or
                    std::is_same_v<Angular<typename F::type>, F>;

} // namespace freq

namespace time
{

template <typename T> struct Seconds : public Unit<Seconds, T> {
    constexpr Seconds(T value) : Unit<Seconds, T>{value} {}
};

template <typename T> struct Samples : public Unit<Samples, T> {
    constexpr explicit Samples(T value) : Unit<Samples, T>{value} {}

    [[nodiscard]] constexpr dsp::freq::NyquistRatio<T> toNyquistRatio() const;
};
} // namespace time

template <typename T>
constexpr time::Samples<T> freq::NyquistRatio<T>::toPeriod() const
{
    return time::Samples{T(2) / this->get()};
}

template <typename T>
constexpr freq::NyquistRatio<T> time::Samples<T>::toNyquistRatio() const
{
    return freq::NyquistRatio{T(2) / this->get()};
}

namespace units
{
constexpr auto operator""_Hz(long double v)
{
    return freq::Hz{static_cast<double>(v)};
}
constexpr auto operator""_Hzf(long double v)
{
    return freq::Hz{static_cast<float>(v)};
}

constexpr auto operator""_rad(long double v)
{
    return freq::Angular{static_cast<double>(v)};
}
constexpr auto operator""_radf(long double v)
{
    return freq::Angular{static_cast<float>(v)};
}

constexpr auto operator""_nyquist(long double v)
{
    return freq::NyquistRatio<double>{static_cast<double>(v)};
}
constexpr auto operator""_nyquistf(long double v)
{
    return freq::NyquistRatio<float>{static_cast<float>(v)};
}

constexpr auto operator""_s(long double v)
{
    return time::Seconds{static_cast<double>(v)};
}
constexpr auto operator""_sf(long double v)
{
    return time::Seconds{static_cast<float>(v)};
}

constexpr auto operator""_ms(long double v)
{
    return time::Seconds{static_cast<double>(v * 0.001l)};
}
constexpr auto operator""_msf(long double v)
{
    return time::Seconds{static_cast<float>(v * 0.001l)};
}

constexpr auto operator""_samples(long double v)
{
    return time::Samples{static_cast<double>(v)};
}
constexpr auto operator""_samplesf(long double v)
{
    return time::Samples{static_cast<float>(v)};
}

using Hz_f32      = dsp::freq::Hz<float>;
using Hz_f64      = dsp::freq::Hz<float>;
using seconds_f32 = dsp::time::Seconds<float>;
using seconds_f64 = dsp::time::Seconds<double>;
} // namespace units

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
