#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "dsp/simd/multi.h"
#include "dsp/simd/simd.h"
#include "dsp/simd/simd_math.h"

using namespace Catch::Matchers;

template <typename T, size_t N> class TestSimd
{
  public:
    static constexpr T kEps = 1e-6;

#define loop(i) for (size_t i = 0; i < N; ++i)

    static void cmp(T x1, T x2)
    {
        if constexpr (!std::is_integral_v<T>) {
            REQUIRE((x1 == x2));
        } else {
            REQUIRE_THAT(x1, WithinAbs(x2, kEps));
        }
    }

    template <int Shift> static void testShift(dsp::simd<T, N> x)
    {
        auto sx = dsp::shift<Shift>(x);

        loop(i)
        {
            auto ii = static_cast<int>(i);
            if ((Shift > 0 && ii < Shift) ||
                (Shift < 0 && ii >= static_cast<int>(N + Shift)))
                cmp(sx[i], 0);
            else {
                cmp(sx[i], x[i - Shift]);
            }
        }
    }

    template <size_t K> static void testPrefix(dsp::simd<T, N> x)
    {
        auto p = dsp::prefix<K>(x);

        T acc[K] = {};
        for (size_t i = 0; i < N; i += K) {
            for (size_t k = 0; k < K; ++k) {
                acc[k] += x[i + k];
                cmp(p[i + k], acc[k]);
            }
        }
    }

    static void run()
    {
        constexpr T kArrayX[] = {4, 2, 39, 19, -15, 17, 23, -100},
                    kArrayY[] = {9, -73, 64, 1, -49, -43, 93, 43};
        auto x                = dsp::simd<T, N>::load(kArrayX);
        auto y                = dsp::simd<T, N>::load(kArrayY);

        SECTION("init")
        {
            T scalar = 1;
            auto z   = dsp::simd<T, N>(scalar);
            loop(i) REQUIRE((z[i] == scalar));
        }
        SECTION("add")
        {
            auto z = x + y;
            loop(i) { cmp(z[i], x[i] + y[i]); }
        }
        SECTION("sub")
        {
            auto z = x - y;
            loop(i) { cmp(z[i], x[i] - y[i]); }
        }
        SECTION("mul")
        {
            if constexpr (!std::is_same_v<T, int64_t>) {
                auto z = x * y;
                loop(i) { cmp(z[i], x[i] * y[i]); }
            }
        }
        if constexpr (!std::is_integral_v<T>) {
            SECTION("div")
            {
                auto z = x / y;
                loop(i) { cmp(z[i], x[i] / y[i]); }
            }
            SECTION("sqrt")
            {
                auto z = sqrt(x);
                loop(i)
                {
                    if (x[i] >= 0) {
                        auto sqrtx = std::sqrt(x[i]);
                        cmp(z[i], sqrtx);
                    }
                }
            }
            SECTION("signbit")
            {
                auto sb = y.signbit();
                loop(i) { cmp(sb[i], std::signbit(y[i])); }
            }
        }
        SECTION("max")
        {
            auto z = max(x, y);
            loop(i) { cmp(z[i], std::max(x[i], y[i])); }
        }
        SECTION("min")
        {
            auto z = min(x, y);
            loop(i) { cmp(z[i], std::min(x[i], y[i])); }
        }
        SECTION("abs")
        {
            auto z = abs(x);
            loop(i) { cmp(z[i], std::abs(x[i])); }
        }
        SECTION("sum")
        {
            auto s1 = sum(x);
            auto s2 = 0;
            loop(i) { s2 += x[i]; }
            cmp(s1, s2);
        }
        SECTION("shift")
        {
            if constexpr (N > 1) {
                testShift<1>(x);
                testShift<-1>(x);
            }
            if constexpr (N > 2) {
                testShift<2>(x);
                testShift<3>(x);
                testShift<-2>(x);
                testShift<-3>(x);
            }
            if constexpr (N > 4) {
                testShift<4>(x);
                testShift<5>(x);
                testShift<6>(x);
                testShift<7>(x);
                testShift<-4>(x);
                testShift<-5>(x);
                testShift<-6>(x);
                testShift<-7>(x);
            }
        }
        SECTION("reduce")
        {
            if constexpr (N > 2) {
                constexpr size_t K = 2;
                auto s1            = dsp::reduce<K>(x);
                T s2[K]{};
                for (size_t k = 0; k < N / K; ++k)
                    for (size_t i = 0; i < K; ++i) s2[i] += x[k * K + i];

                for (size_t i = 0; i < K; ++i) cmp(s1[i], s2[i]);
            }
            if constexpr (N > 4) {
                constexpr size_t K = 4;
                auto s1            = dsp::reduce<K>(x);
                T s2[K]{};
                for (size_t k = 0; k < N / K; ++k)
                    for (size_t i = 0; i < K; ++i) s2[i] += x[k * K + i];

                for (size_t i = 0; i < K; ++i) cmp(s1[i], s2[i]);
            }
        }
        SECTION("flip")
        {
            if constexpr (N > 1) {
                constexpr size_t K = 1;
                auto f             = dsp::flip<K>(x);
                for (size_t i = 0; i < N; ++i) {
                    cmp(f[i], x[i ^ K]);
                }
            }
            if constexpr (N > 2) {
                constexpr size_t K = 2;
                auto f             = dsp::flip<K>(x);
                for (size_t i = 0; i < N; ++i) {
                    cmp(f[i], x[i ^ K]);
                }
            }
            if constexpr (N > 4) {
                constexpr size_t K = 4;
                auto f             = dsp::flip<K>(x);
                for (size_t i = 0; i < N; ++i) {
                    cmp(f[i], x[i ^ K]);
                }
            }
        }
        SECTION("prefix")
        {
            if constexpr (N > 1) {
                testPrefix<1>(x);
            }
            if constexpr (N > 2) {
                testPrefix<2>(x);
            }
            if constexpr (N > 4) {
                testPrefix<4>(x);
            }
        }

        SECTION("push")
        {
            if constexpr (N > 1) {
                auto pval = T(3);
                auto r    = dsp::push(x, pval);
                loop(i)
                {
                    if (i == 0) cmp(r[i], pval);
                    else
                        cmp(r[i], x[i - 1]);
                }
            }
            if constexpr (N > 2) {
                T pvalArr[] = {17, -2};
                auto pval   = dsp::simd<T, 2>::load(pvalArr);
                auto r      = dsp::push(x, pval);
                loop(i)
                {
                    if (i < 2) cmp(r[i], pval[i]);
                    else
                        cmp(r[i], x[i - 2]);
                }
            }
            if constexpr (N > 4) {
                T pvalArr[] = {3, 5, 8, -6};
                auto pval   = dsp::simd<T, 4>::load(pvalArr);
                auto r      = dsp::push(x, pval);
                loop(i)
                {
                    if (i < 4) cmp(r[i], pval[i]);
                    else
                        cmp(r[i], x[i - 4]);
                }
            }
        }

        SECTION("blend")
        {
            auto z = x.cmpgt(y).blend(x, y);
            loop(i) { cmp(z[i], std::max(x[i], y[i])); }
        }

        SECTION("cmpeq")
        {
            auto cmp = x == y;
            loop(i) { REQUIRE((cmp[i] == (x[i] == y[i]))); }
        }
        SECTION("cmpgt")
        {
            auto cmp = x > y;
            loop(i) { REQUIRE((cmp[i] == (x[i] > y[i]))); }
        }
        SECTION("cmpge")
        {
            auto cmp = x >= y;
            loop(i) { REQUIRE((cmp[i] == (x[i] >= y[i]))); }
            cmp = x >= x;
            loop(i) { REQUIRE((cmp[i])); }
        }
        SECTION("cmplt")
        {
            auto cmp = x < y;
            loop(i) { REQUIRE((cmp[i] == (x[i] < y[i]))); }
        }
        SECTION("cmple")
        {
            auto cmp = x <= y;
            loop(i) { REQUIRE((cmp[i] == (x[i] <= y[i]))); }
            cmp = x <= x;
            loop(i) { REQUIRE((cmp[i])); }
        }

        SECTION("any")
        {
            auto any1 = any(x > y);
            bool any2 = false;
            loop(i) { any2 |= x[i] > y[i]; }
            REQUIRE((any1 == any2));

            constexpr T kA[] = {1, 0, 0, 0, 0, 0, 0, 0};
            constexpr T kB[] = {0, 0, 0, 0, 0, 0, 0, 0};
            auto a           = dsp::simd<T, N>::load(kA);
            auto b           = dsp::simd<T, N>::load(kB);
            REQUIRE(any(a > b));
            REQUIRE(any(a >= b));
            REQUIRE(!any(b > a));
        }
        SECTION("all")
        {
            auto all1 = all(x > y);
            bool all2 = true;
            loop(i) { all2 &= x[i] > y[i]; }
            REQUIRE((all1 == all2));

            constexpr T kA[] = {0, 1, 1, 1, 1, 1, 1, 1};
            constexpr T kB[] = {0, 0, 0, 0, 0, 0, 0, 0};
            auto a           = dsp::simd<T, N>::load(kA);
            auto b           = dsp::simd<T, N>::load(kB);
            REQUIRE(!all(a > b));
            REQUIRE(all(a >= b));
            REQUIRE(!all(b > a));
        }

        if constexpr (std::is_integral_v<T>) {
            SECTION("Shift left")
            {
                auto slx = x << 5;
                auto sly = y << 3;
                loop(i)
                {
                    REQUIRE(slx[i] == (x[i] << 5));
                    REQUIRE(sly[i] == (y[i] << 3));
                }
            }
            SECTION("Shift Right")
            {
                auto slx = x >> 8;
                auto sly = y >> 19;
                loop(i)
                {
                    REQUIRE(slx[i] == (x[i] >> 8));
                    REQUIRE(sly[i] == (y[i] >> 19));
                }
            }
        }

        SECTION("store")
        {
            alignas(sizeof(T) * N) T data[16];
            x.store(data);
            loop(i) { REQUIRE(x[i] == data[i]); }
        }
    }
};

TEST_CASE("SIMD float 1", "[dsp][simd]") { TestSimd<float, 1>::run(); }
TEST_CASE("SIMD float 4", "[dsp][simd]") { TestSimd<float, 4>::run(); }
TEST_CASE("SIMD float 2", "[dsp][simd]") { TestSimd<float, 2>::run(); }
#ifdef DSP_SIMD_DOUBLE
TEST_CASE("SIMD double 2", "[dsp][simd]") { TestSimd<double, 2>::run(); }
#endif
TEST_CASE("SIMD int32_t 2", "[dsp][simd]") { TestSimd<int32_t, 2>::run(); }
TEST_CASE("SIMD int32_t 4", "[dsp][simd]") { TestSimd<int32_t, 4>::run(); }
#if DSP_MAX_VEC_SIZE > 16
TEST_CASE("SIMD float 8", "[dsp][simd]") { TestSimd<float, 8>::run(); }
TEST_CASE("SIMD double 4", "[dsp][simd]") { TestSimd<double, 4>::run(); }
TEST_CASE("SIMD int32_t 8", "[dsp][simd]") { TestSimd<int32_t, 8>::run(); }
#endif

TEST_CASE("Convert", "[dsp][simd]")
{
    dsp::mfloat<1> xf1 = {-948232.3239f};
    dsp::mfloat<2> xf2 = {1253.82f, -493.8923f};
    dsp::mfloat<4> xf4 = {8.631f, 2.332f, -39.129f, 20.392f};

    dsp::mint<2> xi2 = {432, -9249};
    dsp::mint<4> xi4 = {19238, -38432, 9383, -492334};

#ifdef DSP_SIMD_DOUBLE
    dsp::mdouble<2> xd2 = {0.5948823, -40.3484834};
#endif

#if DSP_MAX_VEC_SIZE > 16
    dsp::mfloat<8> xf8  = {9394.233f,        0.4943f,        -12.39483f,
                           39334443.f,       -0.0000123929f, 1.239483f,
                           234433.20302023f, 0.233f};
    dsp::mint<8> xi8    = {2938, -3362, 3832, 3, 3234, -393233323, 392, -2833};
    dsp::mdouble<4> xd4 = {0.39848384344334, -393.23343433832, 4.34939494,
                           -393933943.434388383283838};
#endif

    SECTION("as floatx1")
    {
        auto xf2f1 = static_cast<dsp::simd<float, 1>>(xf2);
        REQUIRE((xf2f1[0] == xf2[0]));

        auto xf4f1 = static_cast<dsp::simd<float, 1>>(xf4);
        REQUIRE((xf4f1[0] == xf4[0]));

        auto xi2f1 = static_cast<dsp::simd<float, 1>>(xi2);
        REQUIRE((xi2f1[0] == static_cast<float>(xi2[0])));

        auto xi4f1 = static_cast<dsp::simd<float, 1>>(xi4);
        REQUIRE((xi4f1[0] == static_cast<float>(xi4[0])));

#ifdef DSP_SIMD_DOUBLE
        auto xd2f1 = static_cast<dsp::simd<float, 1>>(xd2);
        REQUIRE((xd2f1[0] == static_cast<float>(xd2[0])));
#endif
    }

    SECTION("as floatx2")
    {
        auto xf1f2 = static_cast<dsp::simd<float, 2>>(xf1);
        REQUIRE((xf1f2[0] == xf1[0]));
        REQUIRE((xf1f2[1] == xf1[0]));

        auto xf4f2 = static_cast<dsp::simd<float, 2>>(xf4);
        REQUIRE((xf4f2[0] == xf4[0]));
        REQUIRE((xf4f2[1] == xf4[1]));

        auto xi2f2 = static_cast<dsp::simd<float, 2>>(xi2);
        REQUIRE((xi2f2[0] == xi2[0]));
        REQUIRE((xi2f2[1] == xi2[1]));

        auto xi4f2 = static_cast<dsp::simd<float, 2>>(xi4);
        REQUIRE((xi4f2[0] == xi4[0]));
        REQUIRE((xi4f2[1] == xi4[1]));

#ifdef DSP_SIMD_DOUBLE
        auto xd2f2 = static_cast<dsp::simd<float, 2>>(xd2);
        REQUIRE((xd2f2[0] == static_cast<float>(xd2[0])));
        REQUIRE((xd2f2[1] == static_cast<float>(xd2[1])));
#endif
    }

    SECTION("as floatx4")
    {
        auto xf1f4 = static_cast<dsp::simd<float, 4>>(xf1);
        REQUIRE((xf1f4[0] == xf1[0]));
        REQUIRE((xf1f4[1] == xf1[0]));
        REQUIRE((xf1f4[2] == xf1[0]));
        REQUIRE((xf1f4[3] == xf1[0]));

        auto xf2f4 = static_cast<dsp::simd<float, 4>>(xf2);
        REQUIRE((xf2f4[0] == xf2[0]));
        REQUIRE((xf2f4[1] == xf2[1]));
        REQUIRE((xf2f4[2] == xf2[0]));
        REQUIRE((xf2f4[3] == xf2[1]));

        auto xi2f4 = static_cast<dsp::simd<float, 4>>(xi2);
        REQUIRE((xi2f4[0] == static_cast<float>(xi2[0])));
        REQUIRE((xi2f4[1] == static_cast<float>(xi2[1])));
        REQUIRE((xi2f4[2] == static_cast<float>(xi2[0])));
        REQUIRE((xi2f4[3] == static_cast<float>(xi2[1])));

        auto xi4f4 = static_cast<dsp::simd<float, 4>>(xi4);
        REQUIRE((xi4f4[0] == static_cast<float>(xi4[0])));
        REQUIRE((xi4f4[1] == static_cast<float>(xi4[1])));
        REQUIRE((xi4f4[2] == static_cast<float>(xi4[2])));
        REQUIRE((xi4f4[3] == static_cast<float>(xi4[3])));

#ifdef DSP_SIMD_DOUBLE
        auto xd2f4 = static_cast<dsp::simd<float, 4>>(xd2);
        REQUIRE((xd2f4[0] == static_cast<float>(xd2[0])));
        REQUIRE((xd2f4[1] == static_cast<float>(xd2[1])));
        REQUIRE((xd2f4[2] == static_cast<float>(xd2[0])));
        REQUIRE((xd2f4[3] == static_cast<float>(xd2[1])));
#endif

#if DSP_MAX_VEC_SIZE > 16
        auto xf8f4 = static_cast<dsp::simd<float, 4>>(xf8);
        REQUIRE((xf8f4[0] == xf8[0]));
        REQUIRE((xf8f4[1] == xf8[1]));
        REQUIRE((xf8f4[2] == xf8[2]));
        REQUIRE((xf8f4[3] == xf8[3]));

        auto xi8f4 = static_cast<dsp::simd<float, 4>>(xi8);
        REQUIRE((xi8f4[0] == static_cast<float>(xi8[0])));
        REQUIRE((xi8f4[1] == static_cast<float>(xi8[1])));
        REQUIRE((xi8f4[2] == static_cast<float>(xi8[2])));
        REQUIRE((xi8f4[3] == static_cast<float>(xi8[3])));

        auto xd4f4 = static_cast<dsp::simd<float, 4>>(xd4);
        REQUIRE((xd4f4[0] == static_cast<float>(xd4[0])));
        REQUIRE((xd4f4[1] == static_cast<float>(xd4[1])));
        REQUIRE((xd4f4[2] == static_cast<float>(xd4[2])));
        REQUIRE((xd4f4[3] == static_cast<float>(xd4[3])));
#endif
    }

    SECTION("as intx2")
    {
        auto xf1i2 = static_cast<dsp::simd<int32_t, 2>>(xf1);
        REQUIRE((xf1i2[0] == static_cast<int32_t>(xf1[0])));
        REQUIRE((xf1i2[1] == static_cast<int32_t>(xf1[0])));

        auto xf2i2 = static_cast<dsp::simd<int32_t, 2>>(xf2);
        REQUIRE((xf2i2[0] == static_cast<int32_t>(xf2[0])));
        REQUIRE((xf2i2[1] == static_cast<int32_t>(xf2[1])));

        auto xi2i2 = static_cast<dsp::simd<int32_t, 2>>(xi2);
        REQUIRE((xi2i2[0] == static_cast<int32_t>(xi2[0])));
        REQUIRE((xi2i2[1] == static_cast<int32_t>(xi2[1])));

        auto xf4i2 = static_cast<dsp::simd<int32_t, 2>>(xf4);
        REQUIRE((xf4i2[0] == static_cast<int32_t>(xf4[0])));
        REQUIRE((xf4i2[1] == static_cast<int32_t>(xf4[1])));

#ifdef DSP_SIMD_DOUBLE
        auto xd2i2 = static_cast<dsp::simd<int32_t, 2>>(xd2);
        REQUIRE((xd2i2[0] == static_cast<int32_t>(xd2[0])));
        REQUIRE((xd2i2[1] == static_cast<int32_t>(xd2[1])));
#endif

#if DSP_MAX_VEC_SIZE > 16
        auto xi8i2 = static_cast<dsp::simd<int32_t, 2>>(xi8);
        REQUIRE((xi8i2[0] == xi8[0]));
        REQUIRE((xi8i2[1] == xi8[1]));

        auto xf8i2 = static_cast<dsp::simd<int32_t, 2>>(xf8);
        REQUIRE((xf8i2[0] == static_cast<int32_t>(xf8[0])));
        REQUIRE((xf8i2[1] == static_cast<int32_t>(xf8[1])));

        auto xd4i2 = static_cast<dsp::simd<int32_t, 2>>(xd4);
        REQUIRE((xd4i2[0] == static_cast<int32_t>(xd4[0])));
        REQUIRE((xd4i2[1] == static_cast<int32_t>(xd4[1])));
#endif
    }

    SECTION("as intx4")
    {
        auto xf1i4 = static_cast<dsp::simd<int32_t, 4>>(xf1);
        REQUIRE((xf1i4[0] == static_cast<int32_t>(xf1[0])));
        REQUIRE((xf1i4[1] == static_cast<int32_t>(xf1[0])));
        REQUIRE((xf1i4[2] == static_cast<int32_t>(xf1[0])));
        REQUIRE((xf1i4[3] == static_cast<int32_t>(xf1[0])));

        auto xf2i4 = static_cast<dsp::simd<int32_t, 4>>(xf2);
        REQUIRE((xf2i4[0] == static_cast<int32_t>(xf2[0])));
        REQUIRE((xf2i4[1] == static_cast<int32_t>(xf2[1])));
        REQUIRE((xf2i4[2] == static_cast<int32_t>(xf2[0])));
        REQUIRE((xf2i4[3] == static_cast<int32_t>(xf2[1])));

        auto xi2i4 = static_cast<dsp::simd<int32_t, 4>>(xi2);
        REQUIRE((xi2i4[0] == static_cast<int32_t>(xi2[0])));
        REQUIRE((xi2i4[1] == static_cast<int32_t>(xi2[1])));
        REQUIRE((xi2i4[2] == static_cast<int32_t>(xi2[0])));
        REQUIRE((xi2i4[3] == static_cast<int32_t>(xi2[1])));

        auto xf4i4 = static_cast<dsp::simd<int32_t, 4>>(xf4);
        REQUIRE((xf4i4[0] == static_cast<int32_t>(xf4[0])));
        REQUIRE((xf4i4[1] == static_cast<int32_t>(xf4[1])));
        REQUIRE((xf4i4[2] == static_cast<int32_t>(xf4[2])));
        REQUIRE((xf4i4[3] == static_cast<int32_t>(xf4[3])));

#ifdef DSP_SIMD_DOUBLE
        auto xd2i4 = static_cast<dsp::simd<int32_t, 4>>(xd2);
        REQUIRE((xd2i4[0] == static_cast<int32_t>(xd2[0])));
        REQUIRE((xd2i4[1] == static_cast<int32_t>(xd2[1])));
        REQUIRE((xd2i4[2] == static_cast<int32_t>(xd2[0])));
        REQUIRE((xd2i4[3] == static_cast<int32_t>(xd2[1])));
#endif

#if DSP_VEC_SIZE > 4
        auto xi8i4 = static_cast<dsp::simd<int32_t, 4>>(xi8);
        REQUIRE((xi8i4[0] == xi8[0]));
        REQUIRE((xi8i4[1] == xi8[1]));
        REQUIRE((xi8i4[2] == xi8[2]));
        REQUIRE((xi8i4[3] == xi8[3]));

        auto xf8i4 = static_cast<dsp::simd<int32_t, 4>>(xf8);
        REQUIRE((xf8i4[0] == static_cast<int32_t>(xf8[0])));
        REQUIRE((xf8i4[1] == static_cast<int32_t>(xf8[1])));
        REQUIRE((xf8i4[2] == static_cast<int32_t>(xf8[2])));
        REQUIRE((xf8i4[3] == static_cast<int32_t>(xf8[3])));

        auto xd4i4 = static_cast<dsp::simd<int32_t, 4>>(xd4);
        REQUIRE((xd4i4[0] == static_cast<int32_t>(xd4[0])));
        REQUIRE((xd4i4[1] == static_cast<int32_t>(xd4[1])));
        REQUIRE((xd4i4[2] == static_cast<int32_t>(xd4[2])));
        REQUIRE((xd4i4[3] == static_cast<int32_t>(xd4[3])));
#endif
    }

#ifdef DSP_SIMD_DOUBLE
    SECTION("as doublex2")
    {
        auto xf1d2 = static_cast<dsp::simd<double, 2>>(xf1);
        REQUIRE((xf1d2[0] == static_cast<double>(xf1[0])));
        REQUIRE((xf1d2[1] == static_cast<double>(xf1[0])));

        auto xf2d2 = static_cast<dsp::simd<double, 2>>(xf2);
        REQUIRE((xf2d2[0] == static_cast<double>(xf2[0])));
        REQUIRE((xf2d2[1] == static_cast<double>(xf2[1])));

        auto xf4d2 = static_cast<dsp::simd<double, 2>>(xf4);
        REQUIRE((xf4d2[0] == static_cast<double>(xf4[0])));
        REQUIRE((xf4d2[1] == static_cast<double>(xf4[1])));

        auto xi2d2 = static_cast<dsp::simd<double, 2>>(xi2);
        REQUIRE((xi2d2[0] == static_cast<double>(xi2[0])));
        REQUIRE((xi2d2[1] == static_cast<double>(xi2[1])));
    }
#endif

#if DSP_VEC_SIZE > 4

    SECTION("as floatx8")
    {
        auto xf1f8 = static_cast<dsp::simd<float, 8>>(xf1);
        REQUIRE((xf1f8[0] == xf1[0]));
        REQUIRE((xf1f8[1] == xf1[0]));
        REQUIRE((xf1f8[2] == xf1[0]));
        REQUIRE((xf1f8[3] == xf1[0]));
        REQUIRE((xf1f8[4] == xf1[0]));
        REQUIRE((xf1f8[5] == xf1[0]));
        REQUIRE((xf1f8[6] == xf1[0]));
        REQUIRE((xf1f8[7] == xf1[0]));

        auto xf2f8 = static_cast<dsp::simd<float, 8>>(xf2);
        REQUIRE((xf2f8[0] == xf2[0]));
        REQUIRE((xf2f8[1] == xf2[1]));
        REQUIRE((xf2f8[2] == xf2[0]));
        REQUIRE((xf2f8[3] == xf2[1]));
        REQUIRE((xf2f8[4] == xf2[0]));
        REQUIRE((xf2f8[5] == xf2[1]));
        REQUIRE((xf2f8[6] == xf2[0]));
        REQUIRE((xf2f8[7] == xf2[1]));

        auto xf4f8 = static_cast<dsp::simd<float, 8>>(xf4);
        REQUIRE((xf4f8[0] == xf4[0]));
        REQUIRE((xf4f8[1] == xf4[1]));
        REQUIRE((xf4f8[2] == xf4[2]));
        REQUIRE((xf4f8[3] == xf4[3]));
        REQUIRE((xf4f8[4] == xf4[0]));
        REQUIRE((xf4f8[5] == xf4[1]));
        REQUIRE((xf4f8[6] == xf4[2]));
        REQUIRE((xf4f8[7] == xf4[3]));

        auto xi2f8 = static_cast<dsp::simd<float, 8>>(xi2);
        REQUIRE((xi2f8[0] == static_cast<float>(xi2[0])));
        REQUIRE((xi2f8[1] == static_cast<float>(xi2[1])));
        REQUIRE((xi2f8[2] == static_cast<float>(xi2[0])));
        REQUIRE((xi2f8[3] == static_cast<float>(xi2[1])));
        REQUIRE((xi2f8[4] == static_cast<float>(xi2[0])));
        REQUIRE((xi2f8[5] == static_cast<float>(xi2[1])));
        REQUIRE((xi2f8[6] == static_cast<float>(xi2[0])));
        REQUIRE((xi2f8[7] == static_cast<float>(xi2[1])));

        auto xi4f8 = static_cast<dsp::simd<float, 8>>(xi4);
        REQUIRE((xi4f8[0] == static_cast<float>(xi4[0])));
        REQUIRE((xi4f8[1] == static_cast<float>(xi4[1])));
        REQUIRE((xi4f8[2] == static_cast<float>(xi4[2])));
        REQUIRE((xi4f8[3] == static_cast<float>(xi4[3])));
        REQUIRE((xi4f8[4] == static_cast<float>(xi4[0])));
        REQUIRE((xi4f8[5] == static_cast<float>(xi4[1])));
        REQUIRE((xi4f8[6] == static_cast<float>(xi4[2])));
        REQUIRE((xi4f8[7] == static_cast<float>(xi4[3])));

        auto xi8f8 = static_cast<dsp::simd<float, 8>>(xi8);
        REQUIRE((xi8f8[0] == static_cast<float>(xi8[0])));
        REQUIRE((xi8f8[1] == static_cast<float>(xi8[1])));
        REQUIRE((xi8f8[2] == static_cast<float>(xi8[2])));
        REQUIRE((xi8f8[3] == static_cast<float>(xi8[3])));
        REQUIRE((xi8f8[4] == static_cast<float>(xi8[4])));
        REQUIRE((xi8f8[5] == static_cast<float>(xi8[5])));
        REQUIRE((xi8f8[6] == static_cast<float>(xi8[6])));
        REQUIRE((xi8f8[7] == static_cast<float>(xi8[7])));

        auto xd2f8 = static_cast<dsp::simd<float, 8>>(xd2);
        REQUIRE((xd2f8[0] == static_cast<float>(xd2[0])));
        REQUIRE((xd2f8[1] == static_cast<float>(xd2[1])));
        REQUIRE((xd2f8[2] == static_cast<float>(xd2[0])));
        REQUIRE((xd2f8[3] == static_cast<float>(xd2[1])));
        REQUIRE((xd2f8[4] == static_cast<float>(xd2[0])));
        REQUIRE((xd2f8[5] == static_cast<float>(xd2[1])));
        REQUIRE((xd2f8[6] == static_cast<float>(xd2[0])));
        REQUIRE((xd2f8[7] == static_cast<float>(xd2[1])));

        auto xd4f8 = static_cast<dsp::simd<float, 8>>(xd4);
        REQUIRE((xd4f8[0] == static_cast<float>(xd4[0])));
        REQUIRE((xd4f8[1] == static_cast<float>(xd4[1])));
        REQUIRE((xd4f8[2] == static_cast<float>(xd4[2])));
        REQUIRE((xd4f8[3] == static_cast<float>(xd4[3])));
        REQUIRE((xd4f8[4] == static_cast<float>(xd4[0])));
        REQUIRE((xd4f8[5] == static_cast<float>(xd4[1])));
        REQUIRE((xd4f8[6] == static_cast<float>(xd4[2])));
        REQUIRE((xd4f8[7] == static_cast<float>(xd4[3])));
    }

    SECTION("as int")
    {
        auto xi2i8 = static_cast<dsp::simd<int32_t, 8>>(xi2);
        REQUIRE((xi2i8[0] == xi2[0]));
        REQUIRE((xi2i8[1] == xi2[1]));
        REQUIRE((xi2i8[2] == xi2[0]));
        REQUIRE((xi2i8[3] == xi2[1]));
        REQUIRE((xi2i8[4] == xi2[0]));
        REQUIRE((xi2i8[5] == xi2[1]));
        REQUIRE((xi2i8[6] == xi2[0]));
        REQUIRE((xi2i8[7] == xi2[1]));

        auto xi4i8 = static_cast<dsp::simd<int32_t, 8>>(xi4);
        REQUIRE((xi4i8[0] == xi4[0]));
        REQUIRE((xi4i8[1] == xi4[1]));
        REQUIRE((xi4i8[2] == xi4[2]));
        REQUIRE((xi4i8[3] == xi4[3]));
        REQUIRE((xi4i8[4] == xi4[0]));
        REQUIRE((xi4i8[5] == xi4[1]));
        REQUIRE((xi4i8[6] == xi4[2]));
        REQUIRE((xi4i8[7] == xi4[3]));

        auto xf4i8 = static_cast<dsp::simd<int32_t, 8>>(xf4);
        REQUIRE((xf4i8[0] == static_cast<int32_t>(xf4[0])));
        REQUIRE((xf4i8[1] == static_cast<int32_t>(xf4[1])));
        REQUIRE((xf4i8[2] == static_cast<int32_t>(xf4[2])));
        REQUIRE((xf4i8[3] == static_cast<int32_t>(xf4[3])));
        REQUIRE((xf4i8[4] == static_cast<int32_t>(xf4[0])));
        REQUIRE((xf4i8[5] == static_cast<int32_t>(xf4[1])));
        REQUIRE((xf4i8[6] == static_cast<int32_t>(xf4[2])));
        REQUIRE((xf4i8[7] == static_cast<int32_t>(xf4[3])));

        auto xf2i8 = static_cast<dsp::simd<int32_t, 8>>(xf2);
        REQUIRE((xf2i8[0] == static_cast<int32_t>(xf2[0])));
        REQUIRE((xf2i8[1] == static_cast<int32_t>(xf2[1])));
        REQUIRE((xf2i8[2] == static_cast<int32_t>(xf2[0])));
        REQUIRE((xf2i8[3] == static_cast<int32_t>(xf2[1])));
        REQUIRE((xf2i8[4] == static_cast<int32_t>(xf2[0])));
        REQUIRE((xf2i8[5] == static_cast<int32_t>(xf2[1])));
        REQUIRE((xf2i8[6] == static_cast<int32_t>(xf2[0])));
        REQUIRE((xf2i8[7] == static_cast<int32_t>(xf2[1])));

        auto xf8i8 = static_cast<dsp::simd<int32_t, 8>>(xf8);
        REQUIRE((xf8i8[0] == static_cast<int32_t>(xf8[0])));
        REQUIRE((xf8i8[1] == static_cast<int32_t>(xf8[1])));
        REQUIRE((xf8i8[2] == static_cast<int32_t>(xf8[2])));
        REQUIRE((xf8i8[3] == static_cast<int32_t>(xf8[3])));
        REQUIRE((xf8i8[4] == static_cast<int32_t>(xf8[4])));
        REQUIRE((xf8i8[5] == static_cast<int32_t>(xf8[5])));
        REQUIRE((xf8i8[6] == static_cast<int32_t>(xf8[6])));
        REQUIRE((xf8i8[7] == static_cast<int32_t>(xf8[7])));

        auto xd2i8 = static_cast<dsp::simd<int32_t, 8>>(xd2);
        REQUIRE((xd2i8[0] == static_cast<int32_t>(xd2[0])));
        REQUIRE((xd2i8[1] == static_cast<int32_t>(xd2[1])));
        REQUIRE((xd2i8[2] == static_cast<int32_t>(xd2[0])));
        REQUIRE((xd2i8[3] == static_cast<int32_t>(xd2[1])));
        REQUIRE((xd2i8[4] == static_cast<int32_t>(xd2[0])));
        REQUIRE((xd2i8[5] == static_cast<int32_t>(xd2[1])));
        REQUIRE((xd2i8[6] == static_cast<int32_t>(xd2[0])));
        REQUIRE((xd2i8[7] == static_cast<int32_t>(xd2[1])));

        auto xd4i8 = static_cast<dsp::simd<int32_t, 8>>(xd4);
        REQUIRE((xd4i8[0] == static_cast<int32_t>(xd4[0])));
        REQUIRE((xd4i8[1] == static_cast<int32_t>(xd4[1])));
        REQUIRE((xd4i8[2] == static_cast<int32_t>(xd4[2])));
        REQUIRE((xd4i8[3] == static_cast<int32_t>(xd4[3])));
        REQUIRE((xd4i8[4] == static_cast<int32_t>(xd4[0])));
        REQUIRE((xd4i8[5] == static_cast<int32_t>(xd4[1])));
        REQUIRE((xd4i8[6] == static_cast<int32_t>(xd4[2])));
        REQUIRE((xd4i8[7] == static_cast<int32_t>(xd4[3])));
    }

    SECTION("as doublex4")
    {
        auto xf1d4 = static_cast<dsp::simd<double, 4>>(xf1);
        REQUIRE((xf1d4[0] == static_cast<double>(xf1[0])));
        REQUIRE((xf1d4[1] == static_cast<double>(xf1[0])));
        REQUIRE((xf1d4[2] == static_cast<double>(xf1[0])));
        REQUIRE((xf1d4[3] == static_cast<double>(xf1[0])));

        auto xf2d4 = static_cast<dsp::simd<double, 4>>(xf2);
        REQUIRE((xf2d4[0] == static_cast<double>(xf2[0])));
        REQUIRE((xf2d4[1] == static_cast<double>(xf2[1])));
        REQUIRE((xf2d4[2] == static_cast<double>(xf2[0])));
        REQUIRE((xf2d4[3] == static_cast<double>(xf2[1])));

        auto xf4d4 = static_cast<dsp::simd<double, 4>>(xf4);
        REQUIRE((xf4d4[0] == static_cast<double>(xf4[0])));
        REQUIRE((xf4d4[1] == static_cast<double>(xf4[1])));
        REQUIRE((xf4d4[2] == static_cast<double>(xf4[2])));
        REQUIRE((xf4d4[3] == static_cast<double>(xf4[3])));

        auto xi2d4 = static_cast<dsp::simd<double, 4>>(xi2);
        REQUIRE((xi2d4[0] == static_cast<double>(xi2[0])));
        REQUIRE((xi2d4[1] == static_cast<double>(xi2[1])));
        REQUIRE((xi2d4[2] == static_cast<double>(xi2[0])));
        REQUIRE((xi2d4[3] == static_cast<double>(xi2[1])));

        auto xi4d4 = static_cast<dsp::simd<double, 4>>(xi4);
        REQUIRE((xi4d4[0] == static_cast<double>(xi4[0])));
        REQUIRE((xi4d4[1] == static_cast<double>(xi4[1])));
        REQUIRE((xi4d4[2] == static_cast<double>(xi4[2])));
        REQUIRE((xi4d4[3] == static_cast<double>(xi4[3])));
    }
#endif
}

TEST_CASE("basic example")
{
    constexpr auto kEps = 1e-9;
    dsp::mfloat<4> x{2.3f, 4.3f, 9.4f, 3.4f};
#ifdef DSP_SIMD_DOUBLE
    dsp::mdouble<2> y{1.4, 9.3};
#else
    dsp::mfloat<2> y{1.4, 9.3};
#endif
    dsp::mint<4> z{3, 4, 5, 6};

    auto res = sqrt(x.load() + y) * z;

    float expected[4];
    expected[0] =
        std::sqrt(x[0] + static_cast<float>(y[0])) * static_cast<float>(z[0]);
    expected[1] =
        std::sqrt(x[1] + static_cast<float>(y[1])) * static_cast<float>(z[1]);
    expected[2] =
        std::sqrt(x[2] + static_cast<float>(y[0])) * static_cast<float>(z[2]);
    expected[3] =
        std::sqrt(x[3] + static_cast<float>(y[1])) * static_cast<float>(z[3]);

    REQUIRE_THAT(res[0], WithinAbs(expected[0], kEps));
    REQUIRE_THAT(res[1], WithinAbs(expected[1], kEps));
    REQUIRE_THAT(res[2], WithinAbs(expected[2], kEps));
    REQUIRE_THAT(res[3], WithinAbs(expected[3], kEps));
}
