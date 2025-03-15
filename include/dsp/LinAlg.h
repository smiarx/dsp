#pragma once

#include "Signal.h"

namespace dsp::linalg
{

/* square matrix */
template <typename F, size_t N>
class Matrix : public std::array<Sample<F, N>, N>
{
  public:
    /* N*N matrix, reprensentation by colum */

    Sample<F, N> mult(Sample<F, N> vec)
    {
        Sample<F, N> y{};
#pragma omp simd
        for (size_t j = 0; j < N; ++j) {
            auto v = vec[j];
            for (size_t i = 0; i < N; ++i) {
                y[i] += v * (*this)[j][i];
            }
        }
        return y;
    }

    Matrix mult(Matrix m)
    {
        Matrix outmat{};
#pragma omp simd
        for (size_t j = 0; j < N; ++j) {
            outmat[j] = mult(Sample<F, N>(m[j]));
        }
        return outmat;
    }
};

template <size_t N> using fMatrix = Matrix<float, N>;

} // namespace dsp::linalg
