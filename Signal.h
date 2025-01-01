#pragma once

#include <array>

namespace dsp
{

template <typename T, int N>
class alignas(N * sizeof(T)) _Signal : public std::array<T, N>
{
    /* parallel _Signal */
  public:
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
