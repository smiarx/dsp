#pragma once

#include "Delay.h"
#include "Lut.h"
#include "Signal.h"
#include <cmath>

namespace dsp
{

namespace kernel
{

template <int A, class Window> class Sinc
{
  public:
    Sinc()                     = delete;
    static constexpr auto Size = A;
    static constexpr auto generate(float x)
    {
        auto xpi = x * M_PIf;
        return sinf(xpi) / (xpi)*Window::generate(x / A);
    }
};

template <int A> class Lanczos
{
  public:
    Lanczos()                  = delete;
    static constexpr auto Size = A;
    static constexpr auto generate(float x)
    {
        auto xpi = x * M_PIf;
        return sinf(xpi) * sinf(xpi / A) / (xpi * xpi) * A;
    }
};
} // namespace kernel

template <int N, class Kernel, int LutSize> class TapKernel : public TapLin<1>
{
  private:
    static constexpr auto A           = Kernel::Size;
    static constexpr auto FilterWidth = A * 2;

    // type used in lookup table
    class KernelType : public std::array<fData<N>, FilterWidth>
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

    /* define the LookUpTable that generate kernels */
    class LutType : public Lut<KernelType, LutSize>
    {
      public:
        LutType()
        {
            // prepare convolution kernel
            Lut<KernelType, LutSize>::fill([](float x) -> auto {
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
                    for (int id = 0; id < FilterWidth; ++id) {
                        auto k     = kernelFromId(id);
                        auto value = Kernel::generate(-k - x);
                        for (int i = 0; i < N; ++i) {
                            kernels[id][i] = value;
                        }
                    }
                }

                return kernels;
            });
        }
    };

  private:
    static LutType lut;

    enum KernelNum {
        First = -A,
        Last  = A - 1,
    };

    static constexpr int kernelFromId(int id) { return First + id; }
    static constexpr int idFromKernel(int k) { return k - First; }

  public:
    template <class Ctxt, class DL> auto read(Ctxt c, const DL &delayline)
    {
        static_assert(Ctxt::VecSize == 1);

        typename Ctxt::Type x = {{0.f}};

        auto idelay = TapNoInterp<1>::id_[0];
        auto fdelay = TapLin<1>::fd_[0];

        auto kernels = lut.read(fdelay);

        constexpr auto VecSize = Ctxt::BaseType::VectorSize;
        for (size_t l = 0; l < FilterWidth - FilterWidth % VecSize;
             l += VecSize) {
            auto points =
                delayline.read(c, idelay - (kernelFromId(l)))[0].toVector();
            inFor(points, k, i) { points[k][i] *= kernels[l + k][i]; }
            inFor(points, k, i) { x[0][i] += points[k][i]; }
        }
        for (int l = FilterWidth - FilterWidth % VecSize; l < FilterWidth;
             ++l) {
            auto point =
                delayline.read(c, idelay - (kernelFromId(l)))[0].toVector();
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
template <int N, class Kernel, int LutSize>
typename TapKernel<N, Kernel, LutSize>::LutType
    TapKernel<N, Kernel, LutSize>::lut;
} // namespace dsp
