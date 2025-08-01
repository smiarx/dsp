#pragma once

#include "Utils.h"
#include "cpu/defines.h"
#include "simd/multi.h"
#include "simd/simd.h"
#include "simd/simd_math.h"
#include <array>

namespace dsp::linalg
{
inline namespace DSP_ARCH_NAMESPACE
{

namespace internal
{

template <class T> struct MatrixInfos {
    static constexpr auto kHeight = 1;
    static constexpr auto kWidth  = 1;
    using type                    = T;
};
template <class M> static constexpr auto kMatrixWidth = MatrixInfos<M>::kWidth;
template <class M>
static constexpr auto kMatrixHeight   = MatrixInfos<M>::kHeight;
template <class M> using MatrixBase_t = typename MatrixInfos<M>::type;

///////////////////////// getSIMD ///////////////////////////////////
template <size_t K, class V> inline auto getSIMD(const V &v)
{
    return v.template getSIMD<K>();
}
template <size_t K> inline auto getSIMD(const float &value) { return value; }
template <size_t K> inline auto getSIMD(const double &value) { return value; }

///////////////////////// Operators ///////////////////////////////////
enum class Op {
    kAdd,
    kSub,
    kMul,
    kDiv,
};

template <Op op, typename X1, typename X2> class VecOp
{
    static_assert((kMatrixHeight<X1> == kMatrixHeight<X2> ||
                   kMatrixHeight<X1> == 1 || kMatrixHeight<X2> == 1) &&
                  kMatrixWidth<X1> == 1 && kMatrixWidth<X2> == 1);

  public:
    VecOp(const X1 &x1, const X2 &x2) : x1_(x1), x2_(x2) {}

    template <size_t K> auto getSIMD() const
    {
        namespace itrn = internal;
        if constexpr (op == Op::kAdd) {
            return itrn::getSIMD<K>(x1_) + itrn::getSIMD<K>(x2_);
        } else if constexpr (op == Op::kSub) {
            return itrn::getSIMD<K>(x1_) - itrn::getSIMD<K>(x2_);
        } else if constexpr (op == Op::kMul) {
            return itrn::getSIMD<K>(x1_) * itrn::getSIMD<K>(x2_);
        } else if constexpr (op == Op::kDiv) {
            return itrn::getSIMD<K>(x1_) / itrn::getSIMD<K>(x2_);
        } else {
            static_assert(false, "Invalid operator");
        }
    }

  private:
    const X1 &x1_;
    const X2 &x2_;
};
template <Op op, typename X1, typename X2>
struct MatrixInfos<VecOp<op, X1, X2>> {
    static constexpr auto kHeight = kMatrixHeight<X1> > kMatrixHeight<X2>
                                        ? kMatrixHeight<X1>
                                        : kMatrixHeight<X2>;
    static constexpr auto kWidth  = 1;
    using type                    = MatrixBase_t<X1>;
};

///////////////////////// AbstractVector ///////////////////////////////////
template <class Vec> class AbstractVector
{
    static_assert(kMatrixWidth<Vec> == 1);

  public:
    AbstractVector(Vec &&vec) : v_(vec) {}

    template <size_t K> inline auto getSIMD() const
    {
        return internal::getSIMD<K>(v_);
    }

    template <class V> auto operator+(const V &other) const
    {
        using namespace internal;
        auto op = VecOp<Op::kAdd, AbstractVector, V>(*this, other);
        return AbstractVector(std::move(op));
    }

    template <class V> auto operator-(const V &other) const
    {
        using namespace internal;
        auto op = VecOp<Op::kSub, AbstractVector, V>(*this, other);
        return AbstractVector(std::move(op));
    }

    template <class V> auto operator*(const V &other) const
    {
        using namespace internal;
        auto op = VecOp<Op::kMul, AbstractVector, V>(*this, other);
        return AbstractVector(std::move(op));
    }

    template <class V> auto operator/(const V &other) const
    {
        using namespace internal;
        auto op = VecOp<Op::kDiv, AbstractVector, V>(*this, other);
        return AbstractVector(std::move(op));
    }

  private:
    const Vec v_;
};

} // namespace internal

///////////////////////// Vector ///////////////////////////////////
template <typename T, size_t N> class Vector
{
    // Vector is store as simd batches,
    // we pad zeros in the front
    // e.g. vector of size 9
    //  {{0 0 a0 a1} {a2 a3 a4} {a5 a6 a7 a8}}
    //

  public:
    static constexpr auto kSIMDSize =
        std::min(DSP_MAX_VEC_SIZE, static_cast<int>(nextPow2(sizeof(T) * N))) /
        sizeof(T);
    static constexpr auto kN     = nextAlignedOffset(N, kSIMDSize);
    static constexpr auto kNsimd = kN / kSIMDSize;

    Vector() = default;
    Vector(const T *data)
    {
        // unaligned
        size_t k = 0;
        if constexpr (kNsimd != kN) {
            for (size_t i = 0; i < kN - N; ++i) {
                data_[k][i] = 0;
            }
            for (size_t i = kN - N; i < kSIMDSize; ++i) {
                data_[k][i] = *data;
                ++data;
            }
            ++k;
        }

        for (; k < kNsimd; ++k) {
            data_[k] = batch<T>::loadu(data);
            data += kSIMDSize;
        }
    }

    T operator[](size_t i)
    {
        i += kN - N;
        auto k = i / kSIMDSize;
        i      = i % kSIMDSize;
        return data_[k][i];
    }

    template <size_t K> [[nodiscard]] auto getSIMD() const
    {
        static_assert(K < kNsimd);
        return data_[K].load();
    }

    template <class V> Vector &operator=(const V &vector)
    {
        store(vector);
        return *this;
    }

    template <class V, size_t K = 0> inline auto store(const V &vector)
    {
        if constexpr (K == kNsimd) return;
        else {
            data_[K].store(vector.template getSIMD<K>());
            store<V, K + 1>(vector);
        }
    }

    template <class V, size_t K = 0> inline auto dot(const V &vector)
    {
        auto res = sum(getSIMD<K>() * (vector.template getSIMD<K>()));
        if constexpr (K + 1 < kNsimd) {
            res += dot<V, K + 1>(vector);
        }
        return res;
    }

    template <class V> auto operator+(const V &other) const
    {
        using namespace internal;
        auto op = VecOp<Op::kAdd, Vector, V>(*this, other);
        return AbstractVector(std::move(op));
    }

    template <class V> auto operator-(const V &other) const
    {
        using namespace internal;
        auto op = VecOp<Op::kSub, Vector, V>(*this, other);
        return AbstractVector(std::move(op));
    }

    template <class V> auto operator*(const V &other) const
    {
        using namespace internal;
        auto op = VecOp<Op::kMul, Vector, V>(*this, other);
        return AbstractVector(std::move(op));
    }

    template <class V> auto operator/(const V &other) const
    {
        using namespace internal;
        auto op = VecOp<Op::kDiv, Vector, V>(*this, other);
        return AbstractVector(std::move(op));
    }

  private:
    std::array<multi<T, kSIMDSize>, kNsimd> data_{};
};

template <typename T, size_t N> struct internal::MatrixInfos<Vector<T, N>> {
    static constexpr auto kHeight = N;
    static constexpr auto kWidth  = 1;
    using type                    = T;
};

///////////////////////// Matrix ///////////////////////////////////
template <typename T, size_t H, size_t W> class Matrix
{
    // Matrix is store as simd batches,
    // we pad zeros in the front rows
    //
  public:
    static constexpr auto kSIMDSize = Vector<T, H>::kSIMDSize;
    static constexpr auto kHsimd    = Vector<T, H>::kNsimd;
    static constexpr auto kNsimd    = kHsimd * W;

    static constexpr auto kSubMatH = kHsimd;
    static constexpr auto kSubMatW = Vector<T, W>::kNsimd;

    static constexpr auto kHOffset = Vector<T, H>::kN - H;
    static constexpr auto kWOffset = Vector<T, W>::kN - W;

    Matrix() = default;
    Matrix(const T *data)
    {
        for (size_t i = 0; i < H; ++i)
            for (size_t j = 0; j < W; ++j) {
                set(i, j, *data);
                ++data;
            }
    }

    void set(size_t i, size_t j, T value)
    {
        i += kHOffset;
        data_[i / kSIMDSize][j][i % kSIMDSize] = value;
    }
    const T &get(size_t i, size_t j) const
    {
        i += kHOffset;
        return data_[i / kSIMDSize][j][i % kSIMDSize];
    }

    template <size_t I, size_t J> auto getSIMD() const
    {
        return data_[I][J].load();
    }

    template <class M, size_t I = 0, size_t J = 0>
    inline void store(const M &matrix)
    {
        if constexpr (J == W) {
            if constexpr (I == kSubMatH) return;
            else {
                store<M, I + 1, 0>(matrix);
            }
        }
        data_[J][I].store(getSIMD<I, J>(matrix));
        store<M, I, J + 1>(matrix);
    }

    template <class V> auto mul(const V &v) const
    {
        return internal::AbstractVector(MatMul(*this, v));
    }

    class MatMul
    {
      public:
        MatMul(const Matrix &mat, const Vector<T, W> &vec) :
            mat_(mat), vec_(vec)
        {
        }

        template <size_t I, size_t K = 0> auto getSIMD() const
        {
            auto vec = vec_.template getSIMD<K>();
            auto res = submatmul<I, K>(vec);

            if constexpr (K + 1 < kSubMatW) {
                res += getSIMD<I, K + 1>();
            }
            return res;
        }

      private:
        template <size_t I, size_t J, size_t K = J == 0 ? kWOffset : 0>
        auto submatmul(simd<T, kSIMDSize> vec) const
        {
            if constexpr (K == kSIMDSize) return T{};
            else {
                auto matsimd =
                    mat_.template getSIMD<I, (J * kSIMDSize - kWOffset + K)>();
                auto res = vec[K] * matsimd;

                if constexpr (K + 1 < kSIMDSize) {
                    res += submatmul<I, J, K + 1>(vec);
                }
                return res;
            };
        };
        const Matrix &mat_;
        const Vector<T, W> &vec_;
    };

  private:
    std::array<std::array<multi<T, kSIMDSize>, W>, kHsimd> data_{};
};

template <typename T, size_t H, size_t W>
struct internal::MatrixInfos<Matrix<T, H, W>> {
    static constexpr auto kHeight = H;
    static constexpr auto kWidth  = W;
    using type                    = T;
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp::linalg
