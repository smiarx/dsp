#pragma once

#include <dsp/cpu/defines.h>
#include <type_traits>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

//=============================================================

template <class T>
concept Numeric = std::is_arithmetic_v<T>;

template <class T>
concept FloatNumeric = std::is_same_v<T, float> || std::is_same_v<T, double>;

//=============================================================

// multi forward declaration
template <typename T, std::size_t N> struct multi;

// ismulti
template <typename T> static constexpr bool kIsMulti = false;
template <typename T, std::size_t N>
static constexpr bool kIsMulti<dsp::multi<T, N>> = true;

template <typename T>
concept IsMulti = kIsMulti<T>;

//=============================================================

// simd forward declaration
template <typename T, std::size_t N> struct simd;
template <typename T, std::size_t N> struct simdmask;

// issimd
template <typename T> constexpr bool kIsSimd                            = false;
template <typename T, std::size_t N> constexpr bool kIsSimd<simd<T, N>> = true;

template <typename T> constexpr bool kIsSimdMask = false;
template <typename T, std::size_t N>
constexpr bool kIsSimdMask<simdmask<T, N>> = true;

template <typename T>
concept IsSimd = kIsSimd<T>;
template <typename T>
concept IsSimdMask = kIsSimdMask<T>;

//=============================================================

template <typename T>
concept MonoSampleType = FloatNumeric<T>;

template <typename T>
concept MultiChannelSampleType =
    (kIsSimd<T> || kIsMulti<T>) && MonoSampleType<typename T::basetype>;

template <typename T>
concept SampleType = MonoSampleType<T> || MultiChannelSampleType<T>;

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
