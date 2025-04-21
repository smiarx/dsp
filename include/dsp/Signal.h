#pragma once

#include "SIMD.h"
#include <array>

namespace dsp
{

/* memory aligned parallel data */
template <typename T = float, size_t N = 1>
class alignas(N * sizeof(T)) Data : public std::array<T, N>
{
  public:
    static_assert(SIMDSIZE % (N * sizeof(T)) == 0,
                  "SIMDSIZE must be divisible by N*sizeof(T)");
    // data type
    using Type = T;

    [[nodiscard]] SIMD_t<T, N> toSIMD() const
    {
        return arrayToSIMD<T, N, true>(this->data());
    }
    void fromSIMD(SIMD_t<T, N> v) { SIMDtoArray<T, N, true>(this->data(), v); }
};

/* N-parallel sample, N should divide VecSize */
template <typename T = float, size_t N = 1> class Sample;

/* L processable N-parallel sample */
template <typename T = float, size_t N = 1, size_t L = 1> class Signal;

template <typename T, size_t N> class Sample : public Data<T, N>
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
    template <bool Vectorize = false> [[nodiscard]] const auto &toSignal() const
    {
        return *reinterpret_cast<
            const std::conditional_t<Vectorize, Vector, Scalar> *>(this);
    }
    auto &toVector() { return toSignal<true>(); }
    [[nodiscard]] const auto &toVector() const { return toSignal<true>(); }
    auto &toScalar() { return toSignal(); }
    [[nodiscard]] const auto &toScalar() const { return toSignal(); }
};

template <typename T, size_t N, size_t L>
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

template <size_t N> using fData = Data<float, N>;
template <size_t N> using dData = Data<double, N>;
template <size_t N> using iData = Data<int, N>;

template <size_t N> using fSample = Sample<float, N>;
template <size_t N> using dSample = Sample<double, N>;
template <size_t N> using iSample = Sample<int, N>;

template <size_t N, size_t L> using fSignal = Signal<float, N, L>;
template <size_t N, size_t L> using dSignal = Signal<double, N, L>;
template <size_t N, size_t L> using iSignal = Signal<int, N, L>;

} // namespace dsp
