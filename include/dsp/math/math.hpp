#pragma once

#include <algorithm>
#include <cmath>

#include <dsp/core/concepts.hpp>
#include <dsp/simd/simd.hpp>

namespace dsp
{

////////////////////// operators ////////////////////////
template <typename T, size_t N>
simd<T, N> always_inline operator+(T s, simd<T, N> x)
{
    return x + s;
}
template <typename T, size_t N>
simd<T, N> always_inline operator*(T s, simd<T, N> x)
{
    return x * s;
}
template <typename T, size_t N>
simd<T, N> always_inline operator-(T s, simd<T, N> x)
{
    return simd<T, N>(s) - x;
}
template <typename T, size_t N>
simd<T, N> always_inline operator/(T s, simd<T, N> x)
{
    return simd<T, N>(s) / x;
}

template <Numeric T> always_inline T get(T x, [[maybe_unused]] size_t i)
{
    return x;
}
template <IsSimd S> always_inline auto get(S x, size_t i) { return x.get(i); }

template <size_t Id, size_t K = 1, typename T> always_inline T getlane(T x)
{
    return x;
}
template <size_t Id, size_t K = 1, IsSimd S> always_inline auto getlane(S x)
{
    return x.template getlane<Id, K>();
}

////////////////////// max ////////////////////////
template <Numeric T> always_inline T max(T x, T y) { return std::max(x, y); }
template <IsSimd S> always_inline S max(S x, S y) { return x.max(y); }

////////////////////// min ////////////////////////
template <Numeric T> always_inline T min(T x, T y) { return std::min(x, y); }
template <IsSimd S> always_inline S min(S x, S y) { return x.min(y); }

////////////////////// abs ////////////////////////
template <Numeric T> always_inline T abs(T x) { return std::abs(x); }
template <IsSimd S> always_inline S abs(S x) { return x.abs(); }

////////////////////// shift ////////////////////////
template <int Shift, Numeric T> always_inline T shift([[maybe_unused]] T x)
{
    return Shift == 0 ? x : T(0);
}
template <int Shift, IsSimd S> always_inline S shift(S x)
{
    return x.template shift<Shift>();
}

////////////////////// sum ////////////////////////
template <Numeric T> always_inline T sum(T x) { return x; }
template <IsSimd T> always_inline auto sum(T x) { return x.sum(); }
////////////////////// product ////////////////////////
template <size_t K = 1, Numeric T> always_inline T product(T x)
{
    static_assert(K == 1);
    return x;
}
template <size_t K = 1, IsSimd S> always_inline auto product(S x)
{
    return x.template product<K>();
}

template <size_t K, Numeric T> always_inline auto reduce(T x)
{
    static_assert(K == 1);
    return x;
}
template <size_t K, IsSimd S> always_inline auto reduce(S x)
{
    return x.template reduce<K>();
}

template <size_t K, IsSimd S> always_inline auto flip(S x)
{
    return x.template flip<K>();
}

template <size_t K = 1, IsSimd S> always_inline auto prefix(S x)
{
    return x.template prefix<K>();
}

/////////////////// duplicate ////////////////////////////
template <size_t K, Numeric T> always_inline auto duplicate(T x)
{
    if constexpr (K == 1) return x;
    else
        return simd<T, K>(x);
}
template <size_t K, IsSimd S> always_inline auto duplicate(S x)
{
    return x.template duplicate<K>();
}

/////////////////// push ////////////////////////////
template <Numeric P, Numeric O>
[[nodiscard]] always_inline auto push(P /*x*/, O y)
{
    return y;
}
template <IsSimd S, typename O>
    requires Numeric<O> or IsSimd<O>
[[nodiscard]] always_inline auto push(S x, O y)
{
    return x.push(y);
}

/////////////////// even/odd ////////////////////////////
template <typename T, size_t N>
[[nodiscard]] always_inline auto even(simd<T, N> x, simd<T, N> y)
{
    return x.even(y);
}
template <typename T, size_t N>
[[nodiscard]] always_inline auto odd(simd<T, N> x, simd<T, N> y)
{
    return x.odd(y);
}

/////////////////// Logic ////////////////////////////

/////////////////// any ////////////////////////////
always_inline bool any(bool x) { return x; }
template <IsSimdMask M> always_inline bool any(M x) { return x.any(); }

/////////////////// all ////////////////////////////
always_inline bool all(bool x) { return x; }
template <IsSimdMask S> always_inline bool all(S x) { return x.all(); }

/////////////////// blend ////////////////////////////
template <Numeric T> always_inline T blend(bool b, T x, T y)
{
    return b ? x : y;
}
template <IsSimdMask M, IsSimd S> always_inline S blend(M m, S x, S y)
{
    return m.blend(x, y);
}

/////////////////// Convert ////////////////////////////

/////////////////// int ////////////////////////////
template <Numeric T> always_inline auto toInt(T x)
{
    return static_cast<int32_t>(x);
}
template <IsSimd S> always_inline auto toInt(S x) { return x.toInt(); }
template <FloatNumeric T> always_inline auto toFloat(int32_t x)
{
    return static_cast<T>(x);
}
template <FloatNumeric T, size_t N>
always_inline auto toFloat(simd<int32_t, N> x)
{
    return x.template toFloat<T>();
}

/////////////////// Math ////////////////////////////

////////////////////// exp ////////////////////////
template <Numeric T> always_inline T exp(T x) { return std::exp(x); }
template <IsSimd S> always_inline S exp(S x) { return x.map(std::exp); }

////////////////////// log ////////////////////////
template <Numeric T> always_inline T log(T x) { return std::log(x); }
template <IsSimd S> always_inline S log(S x) { return x.map(std::log); }

////////////////////// sqrt ////////////////////////
template <Numeric T> always_inline T sqrt(T x) { return std::sqrt(x); }
template <IsSimd S> always_inline S sqrt(S x) { return x.sqrt(); }

////////////////////// signbit ////////////////////////
template <Numeric T> always_inline auto signbit(T x) { return std::signbit(x); }
template <IsSimd S> always_inline auto signbit(S x) { return x.signbit(); }

////////////////////// pow ////////////////////////
template <Numeric T> always_inline T pow(T x, T y) { return std::pow(x, y); }
template <IsSimd S> always_inline S pow(S x, S y) { return x.map(std::pow, y); }
template <Numeric T, size_t N> always_inline simd<T, N> pow(simd<T, N> x, T y)
{
    return x.map(std::pow, simd<T, N>(y));
}
template <Numeric T, size_t N> always_inline simd<T, N> pow(T x, simd<T, N> y)
{
    return simd<T, N>(x).map(std::pow, y);
}

/////////////////// Trig ///////////////////////////

////////////////////// sin ////////////////////////
template <Numeric T> always_inline T sin(T x) { return std::sin(x); }
template <IsSimd S> always_inline S sin(S x) { return x.map(std::sin); }

////////////////////// cos ////////////////////////
template <Numeric T> always_inline T cos(T x) { return std::cos(x); }
template <IsSimd S> always_inline S cos(S x) { return x.map(std::cos); }

////////////////////// tan ////////////////////////
template <Numeric T> always_inline T tan(T x) { return std::tan(x); }
template <IsSimd S> always_inline S tan(S x) { return x.map(std::tan); }

////////////////////// asin ////////////////////////
template <Numeric T> always_inline T asin(T x) { return std::asin(x); }
template <IsSimd S> always_inline S asin(S x) { return x.map(std::asin); }

////////////////////// acos ////////////////////////
template <Numeric T> always_inline T acos(T x) { return std::acos(x); }
template <IsSimd S> always_inline S acos(S x) { return x.map(std::acos); }

////////////////////// atan ////////////////////////
template <Numeric T> always_inline T atan(T x) { return std::atan(x); }
template <IsSimd S> always_inline S atan(S x) { return x.map(std::atan); }

////////////////////// cosh ////////////////////////
template <Numeric T> always_inline T cosh(T x) { return std::cosh(x); }
template <IsSimd S> always_inline S cosh(S x) { return x.map(std::cosh); }

////////////////////// sinh ////////////////////////
template <Numeric T> always_inline T sinh(T x) { return std::sinh(x); }
template <IsSimd S> always_inline S sinh(S x) { return x.map(std::sinh); }

////////////////////// tanh ////////////////////////
template <Numeric T> always_inline T tanh(T x) { return std::tanh(x); }
template <IsSimd S> always_inline S tanh(S x) { return x.map(std::tanh); }

} // namespace dsp
