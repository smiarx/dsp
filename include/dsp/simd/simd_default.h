#pragma once

#include "../cpu/defines.h"
#include <algorithm>
#include <cmath>
#include <cstdint>

#if defined(_MSC_VER)
#define always_inline __forceinline
#define vectorcall    __vectorcall
#else
#define always_inline inline __attribute__((always_inline))
#define vectorcall
#endif

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename> union intFloatUnion;
template <typename> struct floatMask;

template <> union intFloatUnion<float> {
    uint32_t i;
    float f;
};
template <> union intFloatUnion<double> {
    uint64_t i;
    double f;
};

template <> struct floatMask<float> {
    static constexpr intFloatUnion<float> kFull{(uint32_t)-1};
    static constexpr intFloatUnion<float> kSign{0x80000000};
    static constexpr intFloatUnion<float> kNotSign{kFull.i ^ kSign.i};
};

template <> struct floatMask<double> {
    static constexpr intFloatUnion<double> kFull{(uint64_t)-1};
    static constexpr intFloatUnion<double> kSign{0x8000000000000000};
    static constexpr intFloatUnion<double> kNotSign{kFull.i ^ kSign.i};
};

template <typename> struct maskWithType;
template <> struct maskWithType<float> {
    using type = uint32_t;
};
template <> struct maskWithType<double> {
    using type = uint64_t;
};
template <> struct maskWithType<int> {
    using type = uint32_t;
};
template <> struct maskWithType<unsigned int> {
    using type = uint32_t;
};

///////////////////////////////////////////////////////////
// default simd structure using gcc vector extension
template <typename T, size_t N> struct intrin {
#if defined(__GNUC__) || defined(__clang__)
    using basetype                                              = T;
    using type __attribute__((__vector_size__(N * sizeof(T)),
                              __aligned__(N * sizeof(T))))      = T;
    using utype __attribute__((__vector_size__(N * sizeof(T)))) = T;
    using maskbasetype = typename maskWithType<T>::type;
    using masktype __attribute__((__vector_size__(N * sizeof(maskbasetype)),
                                  __aligned__(N * sizeof(maskbasetype)))) =
        maskbasetype;

#define loop(i) for (size_t i = 0; i < N; ++i)

    static constexpr always_inline type init(basetype value)
    {
        type x;
        loop(i) x[i] = value;
        return x;
    }
    static constexpr always_inline type loadu(const basetype *src)
    {
        return *(const utype *)src;
    }
    static constexpr auto load = loadu;
    static void storeu(basetype *x, type v) { *(utype *)x = v; }
    static constexpr auto store = storeu;

    static always_inline type add(type x1, type x2) { return x1 + x2; }
    static always_inline type sub(type x1, type x2) { return x1 - x2; }
    static always_inline type neg(type x) { return -x; }
    static always_inline type mul(type x1, type x2) { return x1 * x2; }
    static always_inline type div(type x1, type x2) { return x1 / x2; }
    static always_inline type sqrt(type x)
    {
        loop(i) x[i] = std::sqrt(x[i]);
        return x;
    }
    static always_inline masktype signbit(type x)
    {
        masktype sb;
        loop(i) sb[i] = std::signbit(x[i]);
        return sb;
    }

    static always_inline masktype bitAnd(masktype x1, masktype x2)
    {
        return x1 & x2;
    }
    static always_inline masktype bitOr(masktype x1, masktype x2)
    {
        return x1 | x2;
    }
    static always_inline masktype bitXor(masktype x1, masktype x2)
    {
        return x1 ^ x2;
    }
    static always_inline masktype bitNot(masktype x) { return ~x; }

    static always_inline type max(type x1, type x2)
    {
        loop(i) x1[i] = std::max(x1[i], x2[i]);
        return x1;
    }
    static always_inline type min(type x1, type x2)
    {
        loop(i) x1[i] = std::min(x1[i], x2[i]);
        return x1;
    }
    static always_inline type abs(type x)
    {
        loop(i) x[i] = std::abs(x[i]);
        return x;
    }
    static always_inline basetype sum(type x)
    {
        basetype s{0};
        loop(i) s += x[i];
        return s;
    }

    static always_inline masktype cmpeq(type x1, type x2)
    {
        masktype mask;
        loop(i) mask[i] = x1[i] == x2[i];
        return mask;
    }
    static always_inline masktype cmpgt(type x1, type x2)
    {
        masktype mask;
        loop(i) mask[i] = x1[i] > x2[i];
        return mask;
    }
    static always_inline masktype cmpge(type x1, type x2)
    {
        masktype mask;
        loop(i) mask[i] = x1[i] >= x2[i];
        return mask;
    }
    static always_inline masktype cmplt(type x1, type x2)
    {
        masktype mask;
        loop(i) mask[i] = x1[i] < x2[i];
        return mask;
    }
    static always_inline masktype cmple(type x1, type x2)
    {
        masktype mask;
        loop(i) mask[i] = x1[i] <= x2[i];
        return mask;
    }
    static always_inline type blend(masktype mask, type x2, type x1)
    {
        loop(i) x1[i] = mask[i] == maskbasetype(0) ? x1[i] : x2[i];
        return x1;
    }

    static always_inline bool any(masktype x)
    {
        bool val = false;
        loop(i) val |= static_cast<bool>(x[i]);
        return val;
    }

    static always_inline bool all(masktype x)
    {
        bool val = true;
        loop(i) val &= static_cast<bool>(x[i]);
        return val;
    }

    template <typename S> static always_inline type convert(S value)
    {
        type retVal;
        using SBase          = decltype(value[0]);
        constexpr auto kSize = sizeof(value) / sizeof(SBase);
        loop(i) { retVal[i] = static_cast<basetype>(value[i % kSize]); }
        return retVal;
    }

#undef loop
#endif
};

/////////////// FLOAT|DOUBLE x 1 //////////////////////////
template <typename T> struct intrin<T, 1> {
    using basetype     = T;
    using type         = T;
    using maskbasetype = bool;
    using masktype     = bool;
    static always_inline type vectorcall init(basetype value) { return value; }
    static always_inline type vectorcall load(const basetype *src)
    {
        return (type)*src;
    }
    static constexpr auto loadu = load;
    static always_inline void vectorcall store(basetype *dest, type value)
    {
        *dest = (basetype)value;
    }
    static constexpr auto storeu = store;
    static always_inline type vectorcall add(type x1, type x2)
    {
        return x1 + x2;
    }
    static always_inline type vectorcall sub(type x1, type x2)
    {
        return x1 - x2;
    }
    static always_inline type vectorcall neg(type x) { return -x; }
    static always_inline type vectorcall mul(type x1, type x2)
    {
        return x1 * x2;
    }
    static always_inline type vectorcall div(type x1, type x2)
    {
        return x1 / x2;
    }
    static always_inline type sqrt(type x) { return std::sqrt(x); }
    static always_inline type signbit(type x) { return std::signbit(x); }
    static always_inline type vectorcall bitAnd(type x1, type x2)
    {
        return x1 & x2;
    }
    static always_inline type vectorcall bitOr(type x1, type x2)
    {
        return x1 | x2;
    }
    static always_inline type vectorcall bitXor(type x1, type x2)
    {
        return x1 ^ x2;
    }
    static always_inline type vectorcall bitNot(type x) { return ~x; }
    static always_inline type vectorcall max(type x1, type x2)
    {
        return std::max(x1, x2);
    }
    static always_inline type vectorcall min(type x1, type x2)
    {
        return std::min(x1, x2);
    }
    static always_inline type vectorcall abs(type x) { return std::abs(x); }

    template <size_t Id, size_t K>
    static always_inline type vectorcall getlane(type x)
    {
        return x;
    }

    static always_inline basetype vectorcall sum(type x) { return x; }

    template <size_t K> static always_inline auto vectorcall duplicate(type x)
    {
        return intrin<T, K>::init(x);
    }

    static always_inline type vectorcall push(type /*x*/, type other)
    {
        return other;
    }

    static always_inline masktype vectorcall cmpeq(type x1, type x2)
    {
        return x1 == x2;
    }
    static always_inline masktype vectorcall cmpgt(type x1, type x2)
    {
        return x1 > x2;
    }
    static always_inline masktype vectorcall cmpge(type x1, type x2)
    {
        return x1 >= x2;
    }
    static always_inline masktype vectorcall cmplt(type x1, type x2)
    {
        return x1 < x2;
    }
    static always_inline masktype vectorcall cmple(type x1, type x2)
    {
        return x1 <= x2;
    }
    static always_inline type vectorcall blend(masktype mask, type x2, type x1)
    {
        return mask == basetype(0) ? x1 : x2;
    }
    static always_inline bool vectorcall any(type mask)
    {
        return mask != basetype(0);
    }
    static always_inline bool vectorcall all(type mask)
    {
        return mask != basetype(0);
    }

    template <typename S> static always_inline type convert(S value)
    {
        return static_cast<type>(value);
    }
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
