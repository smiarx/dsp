#pragma once

#include "Utils.h"
#include "cpu/defines.h"
#include "simd/multi.h"
#include "simd/simd.h"
#include "simd/simd_math.h"
#include <array>
#include <cassert>

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
template <class V>
always_inline auto getSIMD(const V &v, size_t i, size_t j = 0)
{
    return v.getSIMD(i, j);
}
template <>
inline auto getSIMD(const float &value, [[maybe_unused]] size_t i,
                    [[maybe_unused]] size_t j)
{
    return value;
}
template <>
inline auto getSIMD(const double &value, [[maybe_unused]] size_t i,
                    [[maybe_unused]] size_t j)
{
    return value;
}

///////////////////////// Operators ///////////////////////////////////
enum class Op {
    kAdd,
    kSub,
    kMul,
    kDiv,
};

template <Op op, typename X1, typename X2> class MatOp
{
    static_assert((kMatrixHeight<X1> == kMatrixHeight<X2> &&
                   kMatrixWidth<X1> == kMatrixWidth<X2>) ||
                  (kMatrixHeight<X1> == 1 && kMatrixWidth<X1> == 1) ||
                  (kMatrixHeight<X2> == 1 && kMatrixWidth<X2> == 1));

  public:
    MatOp(const X1 &x1, const X2 &x2) : x1_(x1), x2_(x2) {}

    auto getSIMD(size_t i, size_t j = 0) const
    {
        namespace itrn = internal;
        if constexpr (op == Op::kAdd) {
            return itrn::getSIMD(x1_, i, j) + itrn::getSIMD(x2_, i, j);
        } else if constexpr (op == Op::kSub) {
            return itrn::getSIMD(x1_, i, j) - itrn::getSIMD(x2_, i, j);
        } else if constexpr (op == Op::kMul) {
            return itrn::getSIMD(x1_, i, j) * itrn::getSIMD(x2_, i, j);
        } else if constexpr (op == Op::kDiv) {
            return itrn::getSIMD(x1_, i, j) / itrn::getSIMD(x2_, i, j);
        } else {
            static_assert(false, "Invalid operator");
        }
    }

  private:
    const X1 &x1_;
    const X2 &x2_;
};
template <Op op, typename X1, typename X2>
struct MatrixInfos<MatOp<op, X1, X2>> {
    static constexpr auto kHeight = kMatrixHeight<X1> > kMatrixHeight<X2>
                                        ? kMatrixHeight<X1>
                                        : kMatrixHeight<X2>;
    static constexpr auto kWidth  = kMatrixWidth<X1> > kMatrixWidth<X2>
                                        ? kMatrixWidth<X1>
                                        : kMatrixWidth<X2>;
    using type                    = MatrixBase_t<X1>;
};

///////////////////////// AbstractMatrix ///////////////////////////////////
template <class Mat> class AbstractMatrix
{
  public:
    AbstractMatrix(Mat &&mat) : mat_(mat) {}
    AbstractMatrix(const Mat &mat) : mat_(mat) {}

    always_inline auto getSIMD(size_t i, size_t j = 0) const
    {
        return internal::getSIMD(mat_, i, j);
    }

    template <class M> auto operator+(const M &other) const
    {
        using namespace internal;
        auto op = MatOp<Op::kAdd, AbstractMatrix, M>(*this, other);
        return AbstractMatrix<decltype(op)>(std::move(op));
    }

    template <class M> auto operator-(const M &other) const
    {
        using namespace internal;
        auto op = MatOp<Op::kSub, AbstractMatrix, M>(*this, other);
        return AbstractMatrix<decltype(op)>(std::move(op));
    }

    template <class M> auto operator*(const M &other) const
    {
        using namespace internal;
        auto op = MatOp<Op::kMul, AbstractMatrix, M>(*this, other);
        return AbstractMatrix<decltype(op)>(std::move(op));
    }

    template <class M> auto operator/(const M &other) const
    {
        using namespace internal;
        auto op = MatOp<Op::kDiv, AbstractMatrix, M>(*this, other);
        return AbstractMatrix<decltype(op)>(std::move(op));
    }

  private:
    const Mat mat_;
};

} // namespace internal

///////////////////////// Matrix ///////////////////////////////////
template <typename T, size_t H, size_t W> class Matrix
{
    // Matrix is store as simd batches,
    // we pad zeros in the front rows
    //
  public:
    static constexpr auto kSIMDSize =
        std::min(DSP_MAX_VEC_SIZE, static_cast<int>(nextPow2(sizeof(T) * H))) /
        sizeof(T);
    static constexpr auto kH       = nextAlignedOffset(H, kSIMDSize);
    static constexpr auto kSubMatH = kH / kSIMDSize;

    static constexpr auto kW       = nextAlignedOffset(W, kSIMDSize);
    static constexpr auto kSubMatW = kW / kSIMDSize;

    static constexpr auto kNsimd = kSubMatH * W;

    static constexpr auto kHOffset = kH - H;
    static constexpr auto kWOffset = kW - W;

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

    auto getSIMD(size_t i, size_t j) const { return data_[i][j].load(); }
    auto setSIMD(size_t i, size_t j, simd<T, kSIMDSize> s)
    {
        return data_[i][j].store(s);
    }

    template <class M> Matrix &operator=(const M &matrix)
    {
        store(matrix);
        return *this;
    }

    template <class M> inline void store(const M &matrix)
    {
        Matrix out;
        for (size_t i = 0; i < kSubMatH; ++i)
            for (size_t j = 0; j < W; ++j)
                out.setSIMD(i, j, matrix.getSIMD(i, j));
        *this = std::move(out);
    }

    template <class V> auto mul(const V &v) const
    {
        return internal::AbstractMatrix(MatMul(*this, v));
    }

    template <class Mat2> class MatMul
    {
      public:
        MatMul(const Matrix &mat1, const Mat2 &mat2) : mat1_(mat1), mat2_(mat2)
        {
        }

        auto getSIMD(size_t i, size_t j = 0) const
        {
            simd<T, kSIMDSize> out;
            for (size_t k = 0; k < kSubMatW; ++k) {
                auto vec = mat2_.getSIMD(k, j);
                for (size_t l = k == 0 ? kWOffset : 0; l < kSIMDSize; ++l) {
                    auto matsimd =
                        mat1_.getSIMD(i, (k * kSIMDSize - kWOffset + l));
                    auto res = vec[l] * matsimd;
                    if (l == 0 && k == 0) out = res;
                    else
                        out += res;
                }
            }
            return out;
        }

      private:
        const Matrix &mat1_;
        const Mat2 &mat2_;
    };

    template <class M> auto operator+(const M &other) const
    {
        using namespace internal;
        auto op = MatOp<Op::kAdd, Matrix, M>(*this, other);
        return AbstractMatrix<decltype(op)>(std::move(op));
    }

    template <class M> auto operator-(const M &other) const
    {
        using namespace internal;
        auto op = MatOp<Op::kSub, Matrix, M>(*this, other);
        return AbstractMatrix<decltype(op)>(std::move(op));
    }

    template <class M> auto operator*(const M &other) const
    {
        using namespace internal;
        auto op = MatOp<Op::kMul, Matrix, M>(*this, other);
        return AbstractMatrix<decltype(op)>(std::move(op));
    }

    template <class M> auto operator/(const M &other) const
    {
        using namespace internal;
        auto op = MatOp<Op::kDiv, Matrix, M>(*this, other);
        return AbstractMatrix<decltype(op)>(std::move(op));
    }

  private:
    std::array<std::array<multi<T, kSIMDSize>, W>, kSubMatH> data_{};
};

///////////////////////// Vector ///////////////////////////////////
template <typename T, size_t N> class Vector : public Matrix<T, N, 1>
{
    // Vector is store as simd batches,
    // we pad zeros in the front
    // e.g. vector of size 9
    //  {{0 0 a0 a1} {a2 a3 a4} {a5 a6 a7 a8}}
    //

  public:
    using Mat                       = Matrix<T, N, 1>;
    static constexpr auto kSIMDSize = Mat::kSIMDSize;
    static constexpr auto kN        = Mat::kH;
    static constexpr auto kNsimd    = Mat::kSubMatH;

    Vector() = default;
    Vector(const T *data) : Mat(data) {}

    T operator[](size_t i) { return Mat::get(i, 0); }

    [[nodiscard]] auto getSIMD(size_t i, size_t j = 0) const
    {
        return Mat::getSIMD(i, j);
    }

    template <class V> Vector &operator=(const V &vector)
    {
        Mat::store(vector);
        return *this;
    }

    template <class V> inline auto dot(const V &vector)
    {
        T res{};
        for (size_t k = 0; k < kNsimd; ++k)
            res += sum(getSIMD(k) * (vector.getSIMD(k)));
        return res;
    }

    template <class V2> class Outer
    {
        static_assert(internal::kMatrixWidth<V2> == 1);

      public:
        Outer(const Vector &v1, const V2 &v2) : v1_(v1), v2_(v2) {}

        [[nodiscard]] auto getSIMD(size_t i, size_t j) const
        {
            auto v1 = v1_.getSIMD(i);
            j += V2::kHOffset;
            auto v2 = v2_.getSIMD(j / kSIMDSize);

            return v1 * v2[j % kSIMDSize];
        }

      private:
        const Vector &v1_;
        const V2 &v2_;
    };

    template <class V2> auto outer(const V2 &other) const
    {
        return Outer<V2>(*this, other);
    }
};

namespace internal
{
template <typename T, size_t N> struct MatrixInfos<Vector<T, N>> {
    static constexpr auto kHeight = N;
    static constexpr auto kWidth  = 1;
    using type                    = T;
};
} // namespace internal

template <typename T, size_t H, size_t W>
struct internal::MatrixInfos<Matrix<T, H, W>> {
    static constexpr auto kHeight = H;
    static constexpr auto kWidth  = W;
    using type                    = T;
};

} // namespace DSP_ARCH_NAMESPACE
//////////////////////////// Identity ////////////////////////////

template <typename T, size_t N> class Identity
{
  public:
    static constexpr auto kSIMDSize = Matrix<T, N, N>::kSIMDSize;
    static constexpr auto kN        = Matrix<T, N, N>::kH;
    static constexpr auto kSubMatN  = Matrix<T, N, N>::kSubMatH;
    static constexpr auto kNOffset  = Matrix<T, N, N>::kHOffset;
    Identity()                      = default;

    auto getSIMD(size_t i, size_t j) const
    {
        multi<T, kSIMDSize> out{};
        if ((j + kNOffset) / kSIMDSize == i) {
            out[(j + kNOffset) % kSIMDSize] = T{1};
        }
        return out;
    }

    template <class M> auto &mul(const M &mat) { return mat; }

    template <class M> auto operator+(const M &other) const
    {
        using namespace internal;
        auto op = MatOp<Op::kAdd, Identity, M>(*this, other);
        return AbstractMatrix(std::move(op));
    }

    template <class M> auto operator-(const M &other) const
    {
        using namespace internal;
        auto op = MatOp<Op::kSub, Identity, M>(*this, other);
        return AbstractMatrix(std::move(op));
    }

    template <class M> auto operator*(const M &other) const
    {
        using namespace internal;
        auto op = MatOp<Op::kMul, Identity, M>(*this, other);
        return AbstractMatrix(std::move(op));
    }

    template <class M> auto operator/(const M &other) const
    {
        using namespace internal;
        auto op = MatOp<Op::kDiv, Identity, M>(*this, other);
        return AbstractMatrix(std::move(op));
    }
};

} // namespace dsp::linalg
