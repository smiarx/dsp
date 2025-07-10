#pragma once

#include "simd/simd.h"
#include "simd/simd_math.h"
#include <array>

namespace dsp
{

template <typename T, size_t N>
struct alignas(sizeof(T) * N) MultiVal : public std::array<T, N> {

  public:
    static constexpr auto kWidth = N;
    using simdtype               = simd<T, N>;

    always_inline MultiVal()              = default;
    MultiVal(MultiVal &)                  = default;
    MultiVal(const MultiVal &)            = default;
    MultiVal(MultiVal &&)                 = default;
    MultiVal &operator=(const MultiVal &) = default;

    always_inline MultiVal(simdtype value) { value.store((T *)this->data()); }

    template <typename E> constexpr MultiVal(E value) : std::array<T, N>{}
    {
        for (auto &v : *this) v = static_cast<T>(value);
    }
    template <typename... E>
    constexpr MultiVal(E &&...e) :
        std::array<T, N>{{std::forward<const T>(e)...}}
    {
    }

    [[nodiscard]] always_inline static simdtype load(const T *data)
    {
        return simdtype::load(data);
    }

    [[nodiscard]] always_inline simdtype load() const
    {
        return simdtype::load((T *)this->data());
    }

    [[nodiscard]] static always_inline simdtype loadu(const MultiVal *data)
    {
        return simdtype::loadu(data->data());
    }

    always_inline void store(simdtype val) { val.store(this->data()); }

    static always_inline void storeu(MultiVal *dest, simdtype val)
    {
        return val.storeu(dest->data());
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
    always_inline MultiVal &operator*=(MultiVal other)
    {
        store(load() * other.load());
        return *this;
    }

    always_inline MultiVal &operator+=(MultiVal other)
    {
        store(load() + other.load());
        return *this;
    }

    always_inline MultiVal &operator-=(MultiVal other)
    {
        store(load() - other.load());
        return *this;
    }

    always_inline MultiVal &operator/=(MultiVal other)
    {
        store(load() / other.load());
        return *this;
    }

    always_inline auto operator*(MultiVal other) const
    {
        return load() * other.load();
    }

    always_inline auto operator+(MultiVal other) const
    {
        return load() + other.load();
    }

    always_inline auto operator-(MultiVal other) const
    {
        return load() - other.load();
    }

    always_inline auto operator/(MultiVal other) const
    {
        return load() / other.load();
    }

    template <typename T2, size_t N2>
    always_inline MultiVal &operator*=(simd<T2, N2> other)
    {
        store(load() * other);
        return *this;
    }

    template <typename T2, size_t N2>
    always_inline MultiVal &operator+=(simd<T2, N2> other)
    {
        store(load() + other);
        return *this;
    }

    template <typename T2, size_t N2>
    always_inline MultiVal &operator-=(simd<T2, N2> other)
    {
        store(load() - other);
        return *this;
    }

    template <typename T2, size_t N2>
    always_inline MultiVal &operator/=(simd<T2, N2> other)
    {
        store(load() / other);
        return *this;
    }

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

//=============================================================

template <typename T, size_t N>
always_inline auto get(MultiVal<T, N> x, size_t i)
{
    return x[i];
}

//=============================================================

// batch from value
template <typename T> struct Batch {
    using type = MultiVal<T, DSP_MAX_VEC_SIZE / sizeof(T)>;
};

template <typename T, size_t N> struct Batch<MultiVal<T, N>> {
    using type = MultiVal<T, DSP_MAX_VEC_SIZE / sizeof(T)>;
};

template <typename T> using batch = typename Batch<T>::type;

//=============================================================

namespace loadfuncs
{
template <typename T> auto *asBatch(T *ptr)
{
    return reinterpret_cast<batch<T> *>(ptr);
}
template <typename T> const auto *asBatch(const T *ptr)
{
    return reinterpret_cast<const batch<T> *>(ptr);
}

// load and store functions
template <typename T> always_inline T load(T x) { return x; }
template <typename T> always_inline void store(T &dest, T val) { dest = val; }

template <typename T, size_t N>
always_inline auto load(const MultiVal<T, N> &mval)
{
    return mval.load();
}

template <typename T, size_t N>
always_inline auto store(MultiVal<T, N> &dest, simd<T, N> val)
{
    return dest.store(val);
}

template <typename T> auto loadBatch(const T &x)
{
    return batch<T>::loadu(asBatch(&x));
}
template <typename T, typename V> void storeBatch(T &dest, V x)
{
    batch<T>::storeu(asBatch(&dest), x);
}
} // namespace loadfuncs

//=============================================================

// type Width
template <typename T> struct TypeWidth {
    static constexpr auto kWidth = 1;
};
template <typename T, size_t N> struct TypeWidth<simd<T, N>> {
    static constexpr auto kWidth = N;
};
template <typename T, size_t N> struct TypeWidth<MultiVal<T, N>> {
    static constexpr auto kWidth = N;
};
template <typename T> constexpr auto kTypeWidth = TypeWidth<T>::kWidth;

//=============================================================

// int type
template <typename T> struct IntType {
    using type = int;
};
template <typename T, size_t N> struct IntType<MultiVal<T, N>> {
    using type = MultiVal<int, N>;
};
template <typename T, size_t N> struct IntType<simd<T, N>> {
    using type = simd<int, N>;
};
template <typename T> using intType = typename IntType<T>::type;

//=============================================================

// base type
template <typename T> struct BaseType {
    using type = T;
};
template <typename T, size_t N> struct BaseType<simd<T, N>> {
    using type = T;
};
template <typename T, size_t N> struct BaseType<MultiVal<T, N>> {
    using type = T;
};
template <typename T> using baseType = typename BaseType<T>::type;

//=============================================================

// default float, double and int multivals
template <size_t N = (size_t)DSP_MAX_VEC_SIZE / sizeof(float)>
using mfloat = MultiVal<float, N>;

template <size_t N = (size_t)DSP_MAX_VEC_SIZE / sizeof(double)>
using mdouble = MultiVal<double, N>;

template <size_t N = (size_t)DSP_MAX_VEC_SIZE / sizeof(int)>
using mint = MultiVal<int, N>;

} // namespace dsp
