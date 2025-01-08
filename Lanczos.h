#pragma once

#include "Lut.h"
#include "Signal.h"
#include <cmath>

namespace dsp
{
template <int A, int LutSize> class Lanczos
{
    static constexpr auto FilterWidth = A * 2;
    static dsp::Lut<FilterWidth, LutSize> lut_;

    enum Kernel {
        First = -A + 1,
        Last  = A,
    };
    static constexpr int idFromKernel(Kernel k) { return Kernel::Last - k - 1; }

    static constexpr auto kernelFunc(float x)
    {
        auto xpi = x * M_PIf;
        sinf(xpi) * sinf(A * xpi) / (xpi * xpi) * A;
    }

    static void initLut()
    {
        lut_.fill([](float x) -> auto {
            Signal<FilterWidth> kernels = {0.f};
            if (fabs(x) < 1e-7) {
                kernels[A] = 1.f;
            } else {
                for (Kernel k = Kernel::First; k < Kernel::Last; ++k) {
                    auto i = idFromKernel(k);
                    // prepare convolution kernel */
                    kernels[i] = kernelFunc(x - k);
                }
            }
            return kernels;
        });
    }

    template <class Ctxt, class DL>
    auto read(Ctxt c, const DL &delayline, int id, float fd)
    {
        auto kernels = lut_.read(fd);
        auto points =
            delayline.readContiguous<FilterWidth>(c, id - Kernel::Last);

        /* multiply and sum */
        for (int i = 0; i < FilterWidth; ++i) {
            points[i] *= kernels[i];
        }
        auto &x = points[0];
        for (int i = 1; i < FilterWidth; ++i) {
            x += points[i];
        }

        return x;
    }

    template <class Ctxt, class DL>
    auto read(Ctxt c, const DL &delayline, int id, float fd, float scale)
    {
        if (scale < 1.f) {
            return read(c, delayline, id, fd);
        } else {
            // need to scale the kernel by scale
            int Ascaled                      = A * scale;
            decltype(delayline.read(c, 0)) x = {0};
            for (Kernel k = -Ascaled + 1; k < Ascaled; ++k) {
                float kpos  = (k + fd) / scale;
                int kscaled = static_cast<int>(kpos);
                int fkpos   = kpos - kscaled;
                auto kernel = lut_.read(fkpos)[idFromKernel(kscaled)];

                auto point = delayline.read(c, id + k);
                x += point * kernel;
            }
            return x;
        }
    }
};
} // namespace dsp
