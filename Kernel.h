#pragma once

#include "Delay.h"
#include "Lut.h"
#include "Signal.h"
#include <cmath>

namespace dsp
{
template <int N, int A, int LutSize> class TapKernel : public TapLin<1>
{
  private:
    static constexpr auto FilterWidth = A * 2;

    // type used in lookup table
    class KernelType : public std::array<Signal<N>, FilterWidth>
    {
      public:
        KernelType &operator+=(const KernelType &rhs)
        {
            for (int k = 0; k < FilterWidth; ++k) {
                for (int i = 0; i < N; ++i) (*this)[k][i] += rhs[k][i];
            }
            return *this;
        }
        friend auto operator-(KernelType lhs, const KernelType &rhs)
        {
            for (int k = 0; k < FilterWidth; ++k) {
                for (int i = 0; i < N; ++i) lhs[k][i] -= rhs[k][i];
            }
            return lhs;
        }
        friend auto operator*(KernelType lhs, const float rhs)
        {
            for (int k = 0; k < FilterWidth; ++k) {
                for (int i = 0; i < N; ++i) lhs[k][i] *= rhs;
            }
            return lhs;
        }
    };

  public:
    using LutType = Lut<KernelType, LutSize>;

  private:
    static LutType lut;

    enum Kernel {
        First = -A,
        Last  = A - 1,
    };
    static constexpr int kernelFromId(int id) { return First + id; }
    static constexpr int idFromKernel(int k) { return k - First; }

  public:
    static constexpr auto kernelFunc(float x)
    {
        auto xpi = x * M_PIf;
        return sinf(xpi) * sinf(xpi / A) / (xpi * xpi) * A;
    }

  public:
    static void initLut()
    {
        // prepare convolution kernel
        lut.fill([](float x) -> auto {
            KernelType kernels = {{0.f}};
            if (fabs(x) < 1e-7) {
                for (int i = 0; i < N; ++i) {
                    kernels[idFromKernel(0)][i] = 1.f;
                }
            } else if (fabs(x - 1.f) < 1e-7) {
                for (int i = 0; i < N; ++i) {
                    kernels[idFromKernel(-1)][i] = 1.f;
                }
            } else {
                for (int k = Kernel::First; k <= Kernel::Last; ++k) {
                    auto idK   = idFromKernel(k);
                    auto value = kernelFunc(-k - x);
                    for (int i = 0; i < N; ++i) {
                        kernels[idK][i] = value;
                    }
                }
            }
            return kernels;
        });
    }

    template <class Ctxt, class DL> auto read(Ctxt c, const DL &delayline)
    {
        static_assert(Ctxt::VecSize == 1);

        typename Ctxt::Type x = {{0.f}};

        auto id = TapNoInterp<1>::id_[0];
        auto fd = TapLin<1>::fd_[0];

        auto kernels = lut.read(fd);

        constexpr auto VecSize = Ctxt::BaseType::VectorSize;
        for (int l = 0; l < FilterWidth - FilterWidth % VecSize; l += VecSize) {
            auto points =
                delayline.read(c, id - (Kernel::First + l))[0].toVector();
            inFor(points, k, i) { points[k][i] *= kernels[l + k][i]; }
            inFor(points, k, i) { x[0][i] += points[k][i]; }
        }
        for (int l = FilterWidth - FilterWidth % VecSize; l < FilterWidth;
             ++l) {
            auto point =
                delayline.read(c, id - (Kernel::First + l))[0].toVector();
            arrayFor(point[0], i) { x[0][i] += point[0][i] * kernels[l][i]; }
        }

        return x;
    }

    // template <class Ctxt, class DL>
    // auto read(Ctxt c, const DL &delayline, int id, float fd, float scale)
    //{
    //     if (scale < 1.f) {
    //         return read(c, delayline, id, fd);
    //     } else {
    //         // need to scale the kernel by scale
    //         int Ascaled                      = A * scale;
    //         decltype(delayline.read(c, 0)) x = {0};
    //         for (Kernel k = -Ascaled + 1; k < Ascaled; ++k) {
    //             float kpos  = (k + fd) / scale;
    //             int kscaled = static_cast<int>(kpos);
    //             int fkpos   = kpos - kscaled;
    //             auto kernel = lut_.read(fkpos)[idFromKernel(kscaled)];

    //            auto point = delayline.read(c, id + k);
    //            x += point * kernel;
    //        }
    //        return x;
    //    }
    //}
};

// define static variable
template <int N, int A, int LutSize>
typename TapKernel<N, A, LutSize>::LutType TapKernel<N, A, LutSize>::lut;
} // namespace dsp
