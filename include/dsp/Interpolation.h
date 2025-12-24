#pragma once

#include "simd/simd_math.h"
#include <cassert>
#include <cstddef>
#include <tuple>

namespace dsp::interp
{

/*
 *   POLYNOMIALS
 *   ===========
 */
template <typename T, size_t N> class Polynomial : public Polynomial<T, N - 1>
{
  public:
    template <typename... Ts>
    constexpr Polynomial(T p0, Ts... ps) : Polynomial<T, N - 1>(ps...), p_(p0)
    {
    }

    constexpr T operator()(T x) const
    {
        if constexpr (N == 1) {
            return p_;
        } else {
            return p_ + x * Polynomial<T, N - 1>::operator()(x);
        }
    }

    template <size_t D> constexpr T get() const
    {
        // get coefficient of degree D
        if constexpr (D == 0) {
            return p_;
        } else {
            return Polynomial<T, N - 1>::template get<D - 1>();
        }
    }

  private:
    const T p_{};
};

template <typename T> class Polynomial<T, 0>
{
};

template <typename T, typename... Ts>
Polynomial(T, Ts...) -> Polynomial<T, sizeof...(Ts) + 1>;

/*
 * Polynomials extra
 */
template <typename T> struct PolyPoint {
    T x, y;
};
template <typename T>
constexpr auto extrema(const Polynomial<T, 3> &poly) -> PolyPoint<T>
{
    const auto [p0, p1, p2] = poly;
    assert(p2 != 0);

    auto p1div2 = p1 * T(.5);

    auto extrx = -p1div2 / p2;
    auto extry = p0 + p1div2 * extrx;

    return {extrx, extry};
}

template <typename T>
constexpr auto extrema(const Polynomial<T, 4> &poly)
    -> std::tuple<PolyPoint<T>, PolyPoint<T>>
{
    auto [p0, p1, p2, p3] = poly;
    assert(p3 != 0);

    auto p3x3  = 3 * p3;
    auto ip3x3 = T(1) / p3x3;
    auto delta = p2 * p2 - p1 * p3x3;

    assert(delta >= 0);
    auto sqrtDelta = dsp::sqrt(delta);

    auto extr1 = (-p2 - sqrtDelta) * ip3x3;
    auto extr2 = (-p2 + sqrtDelta) * ip3x3;

    return {{extr1, poly(extr1)}, {extr2, poly(extr2)}};
}

/*
 * INTERPOLATION
 * =============
 */

/* parabolic interpolation */
template <typename T> class LagrangeParabolic : public Polynomial<T, 3>
{
    /* Lagrange quadradic interpolation of points
     * at position -1,0,1
     */
  public:
    LagrangeParabolic(T ym1, T y0, T y1) :
        Polynomial<T, 3>{
            y0,
            T(.5) * (y1 - ym1),
            T(.5) * (y1 + ym1) - y0,
        }
    {
    }
};

/* cubic interpolation */
template <typename T> class LagrangeCubic : public Polynomial<T, 4>
{
    /* Lagrange cubic  interpolation of points
     * at position -1,0,1,2
     */
    static constexpr T kOv3 = 0.3333333333333333;
    static constexpr T kOv6 = 0.16666666666666666;

  public:
    LagrangeCubic(T ym1, T y0, T y1, T y2) :
        Polynomial<T, 4>{
            y0,
            -ym1 * kOv3 - T(.5) * y0 + y1 - y2 * kOv6,
            T(.5) * (y1 + ym1) - y0,
            (y2 - ym1) * kOv6 + T(.5) * (y0 - y1),
        }
    {
    }
};

/* hermite interpolation */
template <typename T> class Hermite : public Polynomial<T, 4>
{
    /* Hermite interpolation of points
     * at position -1,0,1,2
     * p(0) = y0, p(1) = y1
     * dp/dx(0) = (y1-ym1)/2, dp/dx(1) = (y2-y0)/2
     */
  public:
    Hermite(T ym1, T y0, T y1, T y2) :
        Polynomial<T, 4>{
            y0,
            (y1 - ym1) * T(.5),
            ym1 - T(2.5) * y0 + T(2) * y1 - y2 * T(.5),
            (y2 - ym1) * T(.5) + T(1.5) * (y0 - y1),
        }
    {
    }
};

} // namespace dsp::interp

// NOLINTBEGIN (readability-identifier-naming)

// structured bindings for polynomial
namespace std
{
template <typename T, size_t N>
struct tuple_size<dsp::interp::Polynomial<T, N>> {
    static constexpr size_t value = N;
};
template <size_t I, typename T, size_t N>
struct tuple_element<I, dsp::interp::Polynomial<T, N>> {
    using type = T;
};
// NOLINTEND (readability-identifier-naming)
} // namespace std
