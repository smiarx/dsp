#pragma once

#include "SIMD.h"
#include <array>

namespace dsp
{

/* memory aligned parallel data */
template <typename T = float, int N = 1>
class alignas(N * sizeof(T)) Data : public std::array<T, N>
{
  public:
    // is using a whole simd register
    static constexpr bool Whole = N % 4 == 0;
    // data type
    using Type = T;

    SIMD_t<T, N> toSIMD() const
    {
        return arrayToSIMD<T, N, Whole>(this->data());
    }
    void fromSIMD(SIMD_t<T, N> v) { SIMDtoArray<T, N, Whole>(this->data(), v); }
};

/* N-parallel sample, N should divide VecSize */
template <typename T = float, int N = 1> class Sample;

/* L processable N-parallel sample */
template <typename T = float, int N = 1, int L = 1> class Signal;

template <typename T, int N> class Sample : public Data<T, N>
{
  public:
    // Ensure that SIMD vectorization is possible
    static_assert(SIMDSIZE % (N * sizeof(T)) == 0,
                  "SIMDSIZE must be divisible by N*sizeof(T)");

    static constexpr auto VectorSize = SIMDSIZE / (N * sizeof(T));
    using Vector                     = Signal<T, N, VectorSize>;
    using Scalar                     = Signal<T, N, 1>;

    template <bool Vectorize = false> auto &toSignal()
    {
        return *reinterpret_cast<
            std::conditional_t<Vectorize, Vector, Scalar> *>(this);
    }
    template <bool Vectorize = false> const auto &toSignal() const
    {
        return *reinterpret_cast<
            const std::conditional_t<Vectorize, Vector, Scalar> *>(this);
    }
    auto &toVector() { return toSignal<true>(); }
    const auto &toVector() const { return toSignal<true>(); }
    auto &toScalar() { return toSignal(); }
    const auto &toScalar() const { return toSignal(); }
};

template <typename T, int N, int L>
class Signal : public std::array<Sample<T, N>, L>
{
  public:
    static constexpr auto Size = N * L;
    static_assert(SIMDSIZE % (Size * sizeof(T)) == 0,
                  "SIMDSIZE must be divisible by signal length");

    SIMD_t<T, Size> toSIMD()
    {
        return arrayToSIMD<T, Size>(this->data()->data());
    }
    void fromSIMD(SIMD_t<T, Size> v)
    {
        SIMDtoArray<T, Size>(this->data()->data(), v);
    }
};

template <int N> using fData = Data<float, N>;
template <int N> using dData = Data<double, N>;
template <int N> using iData = Data<int, N>;

template <int N> using fSample = Sample<float, N>;
template <int N> using dSample = Sample<double, N>;
template <int N> using iSample = Sample<int, N>;

template <int N, int L> using fSignal = Signal<float, N, L>;
template <int N, int L> using dSignal = Signal<double, N, L>;
template <int N, int L> using iSignal = Signal<int, N, L>;

} // namespace dsp
