#pragma once

#include <algorithm>
#include <cmath>

#include "simd.h"

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

template <typename T> always_inline auto get(T x, [[maybe_unused]] size_t i)
{
    return x;
}
template <typename T, size_t N> always_inline auto get(simd<T, N> x, size_t i)
{
    return x.get(i);
}

////////////////////// max ////////////////////////
template <typename T> always_inline T max(T x, T y) { return std::max(x, y); }
template <typename T, size_t N>
always_inline simd<T, N> max(simd<T, N> x, simd<T, N> y)
{
    return x.max(y);
}

////////////////////// min ////////////////////////
template <typename T> always_inline T min(T x, T y) { return std::min(x, y); }
template <typename T, size_t N>
always_inline simd<T, N> min(simd<T, N> x, simd<T, N> y)
{
    return x.min(y);
}

////////////////////// abs ////////////////////////
template <typename T> always_inline T abs(T x) { return std::abs(x); }
template <typename T, size_t N> always_inline simd<T, N> abs(simd<T, N> x)
{
    return x.abs();
}

////////////////////// sum ////////////////////////
template <typename T> always_inline T sum(T x) { return x; }
template <typename T, size_t N> always_inline T sum(simd<T, N> x)
{
    return x.sum();
}
template <size_t K, typename T> always_inline auto reduce(T x)
{
    static_assert(K == 1);
    return x;
}
template <size_t K, typename T, size_t N>
always_inline auto reduce(simd<T, N> x)
{
    return x.template reduce<K>();
}

/////////////////// Logic ////////////////////////////

/////////////////// any ////////////////////////////
template <typename T> always_inline T any(T x) { return x; }
template <typename T, size_t N> always_inline T any(simdmask<T, N> x)
{
    return x.any();
}

/////////////////// all ////////////////////////////
template <typename T> always_inline T all(T x) { return x; }
template <typename T, size_t N> always_inline T all(simdmask<T, N> x)
{
    return x.all();
}

/////////////////// Convert ////////////////////////////

/////////////////// int ////////////////////////////
template <typename T> always_inline auto toInt(T x)
{
    return static_cast<int>(x);
}
template <typename T, size_t N> always_inline auto toInt(simd<T, N> x)
{
    return x.toInt();
}

/////////////////// Math ////////////////////////////

////////////////////// exp ////////////////////////
template <typename T> always_inline T exp(T x) { return std::exp(x); }
template <typename T, size_t N> always_inline simd<T, N> exp(simd<T, N> x)
{
    return x.map(std::exp);
}

////////////////////// log ////////////////////////
template <typename T> always_inline T log(T x) { return std::log(x); }
template <typename T, size_t N> always_inline simd<T, N> log(simd<T, N> x)
{
    return x.map(std::log);
}

////////////////////// sqrt ////////////////////////
template <typename T> always_inline T sqrt(T x) { return std::sqrt(x); }
template <typename T, size_t N> always_inline simd<T, N> sqrt(simd<T, N> x)
{
    return x.sqrt();
}

////////////////////// pow ////////////////////////
template <typename T> always_inline T log(T x, T y) { return std::pow(x, y); }
template <typename T, size_t N>
always_inline simd<T, N> pow(simd<T, N> x, simd<T, N> y)
{
    return x.map(std::pow, y);
}
template <typename T, size_t N> always_inline simd<T, N> pow(simd<T, N> x, T y)
{
    return x.map(std::pow, simd<T, N>::init(y));
}

/////////////////// Trig ///////////////////////////

////////////////////// sin ////////////////////////
template <typename T> always_inline T sin(T x) { return std::sin(x); }
template <typename T, size_t N> always_inline simd<T, N> sin(simd<T, N> x)
{
    return x.map(std::sin);
}

////////////////////// cos ////////////////////////
template <typename T> always_inline T cos(T x) { return std::cos(x); }
template <typename T, size_t N> always_inline simd<T, N> cos(simd<T, N> x)
{
    return x.map(std::cos);
}

////////////////////// tan ////////////////////////
template <typename T> always_inline T tan(T x) { return std::tan(x); }
template <typename T, size_t N> always_inline simd<T, N> tan(simd<T, N> x)
{
    return x.map(std::tan);
}

////////////////////// asin ////////////////////////
template <typename T> always_inline T asin(T x) { return std::asin(x); }
template <typename T, size_t N> always_inline simd<T, N> asin(simd<T, N> x)
{
    return x.map(std::asin);
}

////////////////////// acos ////////////////////////
template <typename T> always_inline T acos(T x) { return std::acos(x); }
template <typename T, size_t N> always_inline simd<T, N> acos(simd<T, N> x)
{
    return x.map(std::acos);
}

////////////////////// atan ////////////////////////
template <typename T> always_inline T atan(T x) { return std::atan(x); }
template <typename T, size_t N> always_inline simd<T, N> atan(simd<T, N> x)
{
    return x.map(std::atan);
}

////////////////////// cosh ////////////////////////
template <typename T> always_inline T cosh(T x) { return std::cosh(x); }
template <typename T, size_t N> always_inline simd<T, N> cosh(simd<T, N> x)
{
    return x.map(std::cosh);
}

////////////////////// sinh ////////////////////////
template <typename T> always_inline T sinh(T x) { return std::sinh(x); }
template <typename T, size_t N> always_inline simd<T, N> sinh(simd<T, N> x)
{
    return x.map(std::sinh);
}

////////////////////// tanh ////////////////////////
template <typename T> always_inline T tanh(T x) { return std::tanh(x); }
template <typename T, size_t N> always_inline simd<T, N> tanh(simd<T, N> x)
{
    return x.map(std::tanh);
}

} // namespace dsp
