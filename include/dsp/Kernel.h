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
    Sinc()                      = delete;
    static constexpr auto kSize = A;
    static constexpr auto generate(float x)
    {
        auto xpi = x * dsp::constants<float>::pi;
        return sinf(xpi) / (xpi)*Window::generate(x / A);
    }
};

template <int A> class Lanczos
{
  public:
    Lanczos()                   = delete;
    static constexpr auto kSize = A;
    static constexpr auto generate(float x)
    {
        auto xpi = x * dsp::constants<float>::pi;
        return sinf(xpi) * sinf(xpi / A) / (xpi * xpi) * A;
    }
};
} // namespace kernel

template <size_t N, class Kernel, size_t LutSize>
class TapKernel : public TapLin<1>
{
  private:
    static constexpr auto kA           = Kernel::kSize;
    static constexpr auto kFilterWidth = kA * 2;

    // type used in lookup table
    class KernelType : public std::array<fData<N>, kFilterWidth>
    {
      public:
        KernelType &operator+=(const KernelType &rhs)
        {
            for (size_t k = 0; k < kFilterWidth; ++k) {
                for (size_t i = 0; i < N; ++i) (*this)[k][i] += rhs[k][i];
            }
            return *this;
        }
        friend auto operator-(KernelType lhs, const KernelType &rhs)
        {
            for (size_t k = 0; k < kFilterWidth; ++k) {
                for (size_t i = 0; i < N; ++i) lhs[k][i] -= rhs[k][i];
            }
            return lhs;
        }
        friend auto operator*(KernelType lhs, const float rhs)
        {
            for (size_t k = 0; k < kFilterWidth; ++k) {
                for (size_t i = 0; i < N; ++i) lhs[k][i] *= rhs;
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
                KernelType kernels = {};
                if (std::fabs(x) < 1e-7) {
                    for (size_t i = 0; i < N; ++i) {
                        kernels[idFromKernel(0)][i] = 1.f;
                    }
                } else if (std::fabs(x - 1.f) < 1e-7) {
                    for (size_t i = 0; i < N; ++i) {
                        kernels[idFromKernel(-1)][i] = 1.f;
                    }
                } else {
                    for (size_t id = 0; id < kFilterWidth; ++id) {
                        auto k     = static_cast<float>(kernelFromId(id));
                        auto value = Kernel::generate(-k - x);
                        for (size_t i = 0; i < N; ++i) {
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
        kFirst = -kA,
        kLast  = kA - 1,
    };

    static constexpr int kernelFromId(size_t id)
    {
        return kFirst + static_cast<int>(id);
    }
    static constexpr size_t idFromKernel(int k)
    {
        return static_cast<size_t>(k - kFirst);
    }

  public:
    template <class Ctxt, class DL> auto read(Ctxt c, const DL &delayline)
    {
        static_assert(Ctxt::kVecSize == 1);

        typename Ctxt::Type x = {};

        auto idelay = TapNoInterp<1>::id_[0];
        auto fdelay = TapLin<1>::fd_[0];

        auto kernels = lut.read(fdelay);

        constexpr auto kVecSize = Ctxt::BaseType::kVectorSize;
        for (size_t l = 0; l < kFilterWidth - kFilterWidth % kVecSize;
             l += kVecSize) {
            auto delay  = idelay - kernelFromId(l);
            auto points = delayline.read(c, delay)[0].toVector();
            inFor(points, k, i) { points[k][i] *= kernels[l + k][i]; }
            inFor(points, k, i) { x[0][i] += points[k][i]; }
        }
        for (size_t l = kFilterWidth - kFilterWidth % kVecSize;
             l < kFilterWidth; ++l) {
            auto delay = idelay - kernelFromId(l);
            auto point = delayline.read(c, delay)[0].toVector();
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
template <size_t N, class Kernel, size_t LutSize>
typename TapKernel<N, Kernel, LutSize>::LutType
    TapKernel<N, Kernel, LutSize>::lut;
} // namespace dsp
