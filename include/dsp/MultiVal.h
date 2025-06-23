#pragma once

#include "simd/simd.h"
#include <array>

namespace dsp
{

template <typename T, size_t N>
struct alignas(sizeof(T) * N) MultiVal : public std::array<T, N> {

  public:
    using simdtype = simd<T, N>;

    always_inline MultiVal() = default;
    always_inline MultiVal(simdtype value) { value.store((T *)this->data()); }
    template <typename... E>
    MultiVal(E &&...e) : std::array<T, N>{{std::forward<T>(e)...}}
    {
    }

    [[nodiscard]] always_inline simdtype load() const
    {
        return simdtype::load((T *)this->data());
    }

    always_inline operator simdtype() const { return load(); }

    template <typename T2, size_t N2>
    always_inline operator simd<T2, N2>() const
    {
        if constexpr (N2 < N) {
            auto reg = simd<T, N2>::load((T *)this->data());
            if constexpr (std::is_same_v<T, T2>) {
                return reg;
            } else {
                return simd<T2, N2>::convert(reg);
            }
        } else {
            return simd<T2, N2>::convert(load());
        }
    }

    // operators
    template <typename T2, size_t N2>
    always_inline auto operator*(simd<T2, N2> other) const
    {
        return other * (*this);
    }

    template <typename T2, size_t N2>
    always_inline auto operator+(simd<T2, N2> other) const
    {
        return other + (*this);
    }

    template <typename T2, size_t N2>
    always_inline auto operator-(simd<T2, N2> other) const
    {
        return static_cast<simd<T2, N2>>(*this) - other;
    }

    template <typename T2, size_t N2>
    always_inline auto operator/(simd<T2, N2> other) const
    {
        return static_cast<simd<T2, N2>>(*this) / other;
    }
};

template <size_t N = (size_t)DSP_VEC_SIZE * 4 / sizeof(float)>
using mfloat = MultiVal<float, N>;

template <size_t N = (size_t)DSP_VEC_SIZE * 4 / sizeof(double)>
using mdouble = MultiVal<double, N>;

template <size_t N = (size_t)DSP_VEC_SIZE * 4 / sizeof(int)>
using mint = MultiVal<int, N>;
} // namespace dsp
