#pragma once

#include "Delay.h"
#include "Lut.h"
#include "Signal.h"
#include <cmath>

namespace dsp
{

namespace kernels
{

template <int A, class Window> class Sinc
{
  public:
    Sinc()                      = delete;
    static constexpr auto kSize = A;
    template <typename F> static constexpr auto generate(F x)
    {
        auto xpi = x * dsp::constants<F>::pi;
        return sin(xpi) / (xpi)*Window::generate(x / A);
    }
};

} // namespace kernels

template <typename T, class Kernel, size_t LutSize>
class TapKernel : public TapLin<baseType<T>>
{
  private:
    using bt                           = baseType<T>;
    static constexpr auto kA           = Kernel::kSize;
    static constexpr auto kFilterWidth = kA * 2;

    // type used in lookup table
    class KernelType : public std::array<T, kFilterWidth>
    {
      public:
        KernelType &operator+=(const KernelType &rhs)
        {
            for (size_t k = 0; k < kFilterWidth; ++k) {
                (*this)[k] += load(rhs[k]);
            }
            return *this;
        }
        friend auto operator-(KernelType lhs, const KernelType &rhs)
        {
            for (size_t k = 0; k < kFilterWidth; ++k) {
                lhs[k] -= load(rhs[k]);
            }
            return lhs;
        }
        friend auto operator*(KernelType lhs, const float rhs)
        {
            for (size_t k = 0; k < kFilterWidth; ++k) {
                lhs[k] = load(lhs[k]) * rhs;
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
                    kernels[idFromKernel(0)] = bt(1);
                } else if (std::fabs(x - 1.f) < 1e-7) {
                    kernels[idFromKernel(-1)] = bt(1);
                } else {
                    for (size_t id = 0; id < kFilterWidth; ++id) {
                        auto k      = static_cast<bt>(kernelFromId(id));
                        auto value  = Kernel::generate(-k - x);
                        kernels[id] = value;
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
        static_assert(!Ctxt::kUseVec);

        decltype(c.getInput()) x{};

        auto idelay = TapNoInterp<bt>::id_;
        auto fdelay = TapLin<bt>::fd_;
        auto delay  = idelay - kernelFromId(0);

        auto kernels = lut.read(fdelay);

        auto l        = 0;
        auto kernCtxt = c.vec();
        kernCtxt.setBlockSize(kFilterWidth);
        kernCtxt.setData(kernels.data());

        CTXTRUN(kernCtxt)
        {
            auto points = delayline.read(kernCtxt, delay);
            points *= kernCtxt.getInput();
            x += reduce<kTypeWidth<T>>(points);
            l += decltype(kernCtxt)::kIncrSize;
        };

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
template <typename T, class Kernel, size_t LutSize>
typename TapKernel<T, Kernel, LutSize>::LutType
    TapKernel<T, Kernel, LutSize>::lut;
} // namespace dsp
