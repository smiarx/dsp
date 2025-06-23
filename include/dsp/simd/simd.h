#pragma once

#include <algorithm>
#include <cstdint>

#include "simd_aarch64.h"
#include "simd_default.h"
#include "simd_x86_64.h"

#ifndef DSP_VEC_SIZE
#define DSP_VEC_SIZE 1
#endif

namespace dsp
{

template <typename T, size_t N> struct simd;

/////////////////// mask ////////////////////////
template <typename T, size_t N> struct simdmask {
    using def          = intrin<T, N>;
    using maskbasetype = typename def::maskbasetype;
    using masktype     = typename def::masktype;

    always_inline simdmask(masktype value) noexcept : value_(value) {}

    always_inline operator masktype() const { return value_; }
    always_inline bool vectorcall any() { return def::any(value_); }

    always_inline simd<T, N> vectorcall blend(simd<T, N> x2, simd<T, N> x1)
    {
        return def::blend(value_, x2, x1);
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

    static always_inline simd load(const T *data) { return def::load(data); }

    always_inline void store vectorcall(basetype *dest) const noexcept
    {
        def::store(dest, value_);
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

    always_inline simd &vectorcall operator&=(simd other) noexcept
    {
        value_ = def::bitAnd(value_, other);
        return *this;
    }

    always_inline simd &vectorcall operator|=(simd other) noexcept
    {
        value_ = def::bitOr(value_, other);
        return *this;
    }

    always_inline simd &vectorcall operator^=(simd other) noexcept
    {
        value_ = def::bitXor(value_, other);
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

    always_inline simd vectorcall operator-() const noexcept
    {
        return def::neg(value_);
    }

    always_inline simd vectorcall operator~() const noexcept
    {
        return def::bitNot(value_);
    }

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

    always_inline simd vectorcall sqrt() const noexcept
    {
        return def::sqrt(value_);
    }

    always_inline basetype vectorcall sum() const noexcept
    {
        return def::sum(value_);
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

  private:
    type value_;
};

} // namespace dsp
