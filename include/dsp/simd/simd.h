#pragma once

#include <algorithm>
#include <cstdint>

#include "../cpu/defines.h"
#include "simd_arm.h"
#include "simd_default.h"
#include "simd_x86_64.h"

#ifndef DSP_MAX_VEC_SIZE
#define DSP_MAX_VEC_SIZE 16
#endif

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T, size_t N> struct simd;

/////////////////// mask ////////////////////////
// mask tool
template <typename T> struct MaskBoolType;
template <> struct MaskBoolType<bool> {
    using type = bool;
};
template <> struct MaskBoolType<int32_t> {
    using type = uint32_t;
};
template <> struct MaskBoolType<uint32_t> {
    using type = uint32_t;
};
template <> struct MaskBoolType<uint64_t> {
    using type = uint64_t;
};
template <> struct MaskBoolType<float> {
    using type = uint32_t;
};
template <> struct MaskBoolType<double> {
    using type = int64_t;
};

template <typename T, size_t N> struct simdmask {
    using def          = intrin<T, N>;
    using maskbasetype = typename def::maskbasetype;
    using masktype     = typename def::masktype;

    always_inline simdmask(masktype value) noexcept : value_(value) {}
    always_inline simdmask() noexcept : value_{} {}

    always_inline operator masktype() const { return value_; }

    always_inline bool vectorcall any() { return def::any(value_); }
    always_inline bool vectorcall all() { return def::all(value_); }

    always_inline bool vectorcall operator[](size_t index) const noexcept
    {
        return get(index);
    }
    union simdunion {
        masktype simd;
        alignas(sizeof(maskbasetype) * N)
            typename MaskBoolType<maskbasetype>::type scalar[N];
    };

    [[nodiscard]] always_inline bool vectorcall get(size_t index) const noexcept
    {
        simdunion u{value_};
        return u.scalar[index];
    }

    always_inline simdmask vectorcall operator&&(simdmask other)
    {
        return def::bitAnd(value_, other);
    }
    always_inline simdmask vectorcall operator||(simdmask other)
    {
        return def::bitOr(value_, other);
    }
    always_inline simdmask vectorcall operator^(simdmask other)
    {
        return def::bitXor(value_, other);
    }

    always_inline simdmask vectorcall operator~() const noexcept
    {
        return def::bitNot(value_);
    }

    always_inline simd<T, N> vectorcall blend(simd<T, N> x2, simd<T, N> x1)
    {
        return def::blend(value_, x2, x1);
    }
    template <typename T2>
    always_inline simd<T2, N> vectorcall blend(simd<T2, N> x2, simd<T2, N> x1)
    {
        union {
            masktype t1;
            simdmask<T2, N> t2;
        } value{value_};
        return intrin<T2, N>::blend(value.t2, x2, x1);
    }

  private:
    masktype value_;
};

////////////////// simd /////////////////////////

// simd type
template <typename T, size_t N> struct simd {
    using def      = intrin<T, N>;
    using type     = typename def::type;
    using basetype = typename def::basetype;
    using mask     = simdmask<T, N>;

    static constexpr auto kWidth = N;

    union simdunion {
        type simd;
        alignas(sizeof(basetype) * N) basetype scalar[N];
    };

    always_inline simd() noexcept : value_(def::init(basetype(0))) {}
    always_inline simd(type value) noexcept : value_(value) {}
    // always_inline simd(std::initializer_list<T> values) noexcept
    //{
    //     value_ = def::loadu(data(values));
    // }

    // scalar constructor (only if type isn't already scalar)
    template <size_t N2 = N, typename = typename std::enable_if_t<(N2 > 1)>>
    always_inline simd(basetype scalar) noexcept : value_(def::init(scalar))
    {
    }

    always_inline operator type() const { return value_; }

    template <typename T2, size_t N2>
    static always_inline simd convert(simd<T2, N2> value)
    {
        return def::convert(value);
    }
    static always_inline simd convert(simd<T, N> value) { return value; }

    [[nodiscard]] static always_inline simd load(const T *data)
    {
        return def::load(data);
    }
    [[nodiscard]] static always_inline simd loadu(const T *data)
    {
        return def::loadu(data);
    }

    always_inline void store vectorcall(basetype *dest) const noexcept
    {
        def::store(dest, value_);
    }

    always_inline void storeu vectorcall(basetype *dest) const noexcept
    {
        def::storeu(dest, value_);
    }

    always_inline basetype vectorcall get(size_t index) const noexcept
    {
        simdunion u{value_};
        return u.scalar[index];
    }

    always_inline basetype vectorcall operator[](size_t index) const noexcept
    {
        return get(index);
    }

    always_inline void vectorcall set(size_t index, basetype scalar) noexcept
    {
        simdunion u{value_};
        u.scalar[index] = scalar;
        value_          = u.simd;
    }

    // increment operators

    always_inline simd &vectorcall operator+=(simd other) noexcept
    {
        value_ = def::add(value_, other);
        return *this;
    }

    always_inline simd &vectorcall operator-=(simd other) noexcept
    {
        value_ = def::sub(value_, other);
        return *this;
    }

    always_inline simd &vectorcall operator*=(simd other) noexcept
    {
        value_ = def::mul(value_, other);
        return *this;
    }

    always_inline simd &vectorcall operator/=(simd other) noexcept
    {
        value_ = def::div(value_, other);
        return *this;
    }

    always_inline simd &vectorcall operator+=(basetype scalar) noexcept
    {
        value_ = def::add(value_, def::init(scalar));
        return *this;
    }

    always_inline simd &vectorcall operator-=(basetype scalar) noexcept
    {
        value_ = def::sub(value_, def::init(scalar));
        return *this;
    }

    always_inline simd &vectorcall operator*=(basetype scalar) noexcept
    {
        value_ = def::mul(value_, def::init(scalar));
        return *this;
    }

    // operators

    always_inline simd vectorcall operator+(simd other) const noexcept
    {
        return def::add(value_, other);
    }

    always_inline simd vectorcall operator-(simd other) const noexcept
    {
        return def::sub(value_, other);
    }

    always_inline simd vectorcall operator*(simd other) const noexcept
    {
        return def::mul(value_, other);
    }

    always_inline simd vectorcall operator/(simd other) const noexcept
    {
        return def::div(value_, other);
    }

    always_inline simd vectorcall operator&(simd other) const noexcept
    {
        return def::bitAnd(value_, other);
    }

    always_inline simd vectorcall operator|(simd other) const noexcept
    {
        return def::bitOr(value_, other);
    }

    always_inline simd vectorcall operator^(simd other) const noexcept
    {
        return def::bitXor(value_, other);
    }

    always_inline simd vectorcall operator<<(int shift) const noexcept
    {
        return def::bitShiftLeft(value_, shift);
    }
    always_inline simd vectorcall operator>>(int shift) const noexcept
    {
        return def::bitShiftRight(value_, shift);
    }

    always_inline simd vectorcall operator-() const noexcept
    {
        return def::neg(value_);
    }

    always_inline mask vectorcall operator<(simd other) const noexcept
    {
        return def::cmplt(value_, other);
    }

    always_inline mask vectorcall operator>(simd other) const noexcept
    {
        return def::cmpgt(value_, other);
    }

    always_inline mask vectorcall operator<=(simd other) const noexcept
    {
        return def::cmple(value_, other);
    }

    always_inline mask vectorcall operator>=(simd other) const noexcept
    {
        return def::cmpge(value_, other);
    }

    always_inline mask vectorcall operator==(simd other) const noexcept
    {
        return def::cmpeq(value_, other);
    }

    // scalar operators

    always_inline simd vectorcall operator+(basetype scalar) const noexcept
    {
        return def::add(value_, def::init(scalar));
    }

    always_inline simd vectorcall operator-(basetype scalar) const noexcept
    {
        return def::sub(value_, def::init(scalar));
    }

    always_inline simd vectorcall operator*(basetype scalar) const noexcept
    {
        return def::mul(value_, def::init(scalar));
    }

    always_inline simd vectorcall operator/(basetype scalar) const noexcept
    {
        return def::div(value_, def::init(scalar));
    }

    always_inline simd vectorcall operator&(simd other)
    {
        return def::bitAnd(value_, other);
    }
    always_inline simd vectorcall operator|(simd other)
    {
        return def::bitOr(value_, other);
    }
    always_inline simd vectorcall operator^(simd other)
    {
        return def::bitXor(value_, other);
    }

    always_inline simd vectorcall sqrt() const noexcept
    {
        return def::sqrt(value_);
    }

    always_inline mask vectorcall signbit() const noexcept
    {
        return def::signbit(value_);
    }

    always_inline basetype vectorcall sum() const noexcept
    {
        return def::sum(value_);
    }

    template <size_t K>
    always_inline std::conditional_t<K == 1, basetype, simd<T, K>>
        vectorcall reduce() const noexcept
    {
        static_assert(K >= 1 && K <= N);
        if constexpr (K == 1) return sum();
        else if constexpr (K == N) {
            return *this;
        } else if constexpr (K == 2) {
            return def::reduce2(value_);
        } else if constexpr (K == 4) {
            return def::reduce4(value_);
        }
    }
    template <size_t K> always_inline simd vectorcall flip() const noexcept
    {
        static_assert(K >= 1 && K <= N);
        if constexpr (K == 1) {
            return def::flip1(value_);
        } else if constexpr (K == 2) {
            return def::flip2(value_);
        } else if constexpr (K == 4) {
            return def::flip4(value_);
        }
    }

    always_inline simd vectorcall max(simd other)
    {
        return def::max(value_, other);
    }

    always_inline simd vectorcall min(simd other)
    {
        return def::min(value_, other);
    }

    always_inline simd vectorcall abs() { return def::abs(value_); }

    always_inline mask vectorcall cmpeq(simd other)
    {
        return def::cmpeq(value_, other);
    }

    always_inline mask vectorcall cmpgt(simd other)
    {
        return def::cmpgt(value_, other);
    }

    always_inline simd map(T (*func)(T))
    {
        simdunion u{value_};
        for (size_t i = 0; i < N; ++i) {
            u.scalar[i] = func(u.scalar[i]);
        }
        return u.simd;
    }

    always_inline simd map(T (*func)(T, T), simd other)
    {
        simdunion u1{value_}, u2{other};
        for (size_t i = 0; i < N; ++i) {
            u1.scalar[i] = func(u1.scalar[i], u2.scalar[i]);
        }
        return u1.simd;
    }

    always_inline auto vectorcall toInt()
    {
        return simd<int32_t, N>::convert(*this);
    }
    template <typename T2> always_inline auto vectorcall toFloat()
    {
        return simd<T2, N>::convert(*this);
    }

  private:
    type value_;
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
