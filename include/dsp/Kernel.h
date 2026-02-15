#pragma once

#include "Delay.h"
#include "Lut.h"
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

    static constexpr auto kMaxScale = 3;

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
        template <typename Float>
        friend auto operator*(KernelType lhs, const Float rhs)
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
                    T sum{};
                    for (size_t id = 0; id < kFilterWidth; ++id) {
                        auto k      = static_cast<bt>(kernelFromId(id));
                        auto value  = Kernel::generate(-k - x);
                        kernels[id] = value;
                        sum += value;
                    }
                    for (size_t id = 0; id < kFilterWidth; ++id) {
                        kernels[id] /= sum;
                    }
                }

                return kernels;
            });
        }
    };

  private:
    LutType lut_;

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

        auto kernels = lut_.read(fdelay);

        auto kernCtxt = c.vec();
        kernCtxt.setBlockSize(kFilterWidth);
        kernCtxt.setData(kernels.data());

        CTXTRUN(kernCtxt)
        {
            auto points = delayline.read(kernCtxt, delay);
            points *= kernCtxt.getInput();
            x += reduce<kTypeWidth<T>>(points);
        };

        return x;
    }

    template <class Ctxt, class DL>
    auto read(Ctxt c, const DL &delayline, bt scale)
    {
        /* bandlimited read delayline with cutoff defined by
         * `NyquistFreq/scale`.
         * to be used if delayline is read at a rate of `scale`
         * much slower than simple read
         */
        static_assert(!Ctxt::kUseVec);

        if (scale <= bt(1)) return read(c, delayline);

        scale         = dsp::min(scale, bt(kMaxScale));
        auto invscale = bt(1) / scale;

        decltype(c.getInput()) x{};

        auto idelay = TapNoInterp<bt>::id_;
        auto fdelay = TapLin<bt>::fd_;

        std::array<T, kFilterWidth * kMaxScale> kernels;

        auto delaywidth = static_cast<int>(fdelay + kA * scale);
        auto pos        = (fdelay - bt(delaywidth)) * invscale;
        size_t i        = 0;
        while (pos < kA * bt(0.999)) {
            int kernel = static_cast<int>(std::floor(pos));
            auto posf  = pos - bt(kernel);
            kernels[i] = lut_.read(posf)[idFromKernel(kernel)];
            ++i;
            pos += invscale;
        }

        auto delay  = idelay + delaywidth;
        auto length = i;

        auto kernCtxt = c.vec();
        kernCtxt.setBlockSize(static_cast<int>(length));
        kernCtxt.setData(kernels.data());

        T kernsum{};
        CTXTRUN(kernCtxt)
        {
            auto points = delayline.read(kernCtxt, delay);
            auto kerns  = kernCtxt.getInput();
            points *= kerns;
            x += reduce<kTypeWidth<T>>(points);
            kernsum += reduce<kTypeWidth<T>>(kerns);
        };

        return x / kernsum;
    }
};
} // namespace dsp
