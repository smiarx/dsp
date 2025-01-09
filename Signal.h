#pragma once

#include <array>

#ifndef SIMDSIZE
#if defined(__AVX2__)
#define SIMDSIZE 8
#elif defined(__SSE__)
#define SIMDSIZE 4
#else
#define SIMDSIZE 1
#endif
#endif

namespace dsp
{

template <typename T, int N>
class alignas(N * sizeof(T)) _Signal : public std::array<T, N>
{
    /* parallel Signal */
  public:
    static constexpr auto VectorSize = SIMDSIZE / N;

    /* when an array of signals is given, use this to get as many signals as
     * possible in a vector */
    using Vector =
        std::array<std::enable_if_t<SIMDSIZE % N == 0, _Signal>, VectorSize>;
    using Scalar = std::array<_Signal, 1>;

    template <class SV> auto &to() { return *reinterpret_cast<SV *>(this); }
    template <class SV> const auto &to() const
    {
        return *reinterpret_cast<const SV *>(this);
    }
    auto &toVector() { return to<Vector>(); }
    const auto &toVector() const { return to<Vector>(); }

    auto &toScalar() { return to<Scalar>(); }
    const auto &toScalar() const { return to<Scalar>(); }

    _Signal &operator=(const T &rhs)
    {
        for (int i = 0; i < N; ++i) (*this)[i] = rhs;
        return *this;
    }

    _Signal &operator+=(const _Signal &rhs)
    {
        for (int i = 0; i < N; ++i) (*this)[i] += rhs[i];
        return *this;
    }

    friend _Signal operator+(_Signal lhs, const _Signal &rhs)
    {
        for (int i = 0; i < N; ++i) lhs[i] += rhs[i];
        return lhs;
    }

    _Signal &operator-=(const _Signal &rhs)
    {
        for (int i = 0; i < N; ++i) this->x[i] -= rhs[i];
        return *this;
    }

    friend _Signal operator-(_Signal lhs, const _Signal &rhs)
    {
        for (int i = 0; i < N; ++i) lhs[i] -= rhs[i];
        return lhs;
    }

    _Signal &operator*=(const _Signal &rhs)
    {
        for (int i = 0; i < N; ++i) (*this)[i] *= rhs[i];
        return *this;
    }

    friend _Signal operator*(_Signal lhs, const _Signal &rhs)
    {
        for (int i = 0; i < N; ++i) lhs[i] *= rhs[i];
        return lhs;
    }

    friend _Signal operator*(_Signal lhs, const T &rhs)
    {
        for (int i = 0; i < N; ++i) lhs[i] *= rhs;
        return lhs;
    }

    friend _Signal operator*(const T &lhs, _Signal rhs)
    {
        for (int i = 0; i < N; ++i) rhs[i] *= lhs;
        return rhs;
    }

    friend _Signal operator/(_Signal lhs, const T &rhs)
    {
        for (int i = 0; i < N; ++i) lhs[i] /= rhs;
        return lhs;
    }

    friend _Signal operator/(const T &lhs, _Signal rhs)
    {
        for (int i = 0; i < N; ++i) rhs[i] /= lhs;
        return rhs;
    }
};

template <int N> using Signal = _Signal<float, N>;

template <int N> using iSignal = _Signal<int, N>;

} // namespace dsp
