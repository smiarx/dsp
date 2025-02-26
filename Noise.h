#pragma once

#include "Signal.h"
#include <cstdlib>
#include <limits>

namespace dsp
{

template <int N> class iNoise
{
    static constexpr int generator[] = {1736462733, 1487623987, 848278349,
                                        1987647829, 1837654627, 963782763,
                                        1492758273, 1746273223};

  public:
    iNoise()
    {
        for (int i = 0; i < N; ++i) {
            seed_[i] = rand();
        }
    }

    auto process()
    {
        for (int i = 0; i < N; ++i) {
            state_[i] = seed_[i] + state_[i] * generator[i];
        }
        return state_;
    }

  private:
    iData<N> seed_;
    iData<N> state_{0};
};

template <int N> class Noise : public iNoise<N>
{
  public:
    auto process()
    {
        using type  = typename fData<N>::Type;
        using itype = typename iData<N>::Type;
        fData<N> y;

        auto iy = iNoise<N>::process();
        for (int i = 0; i < N; ++i) {
            y[i] = static_cast<type>(iy[i]) /
                   static_cast<type>(std::numeric_limits<itype>::max());
        }
        return y;
    }
};

} // namespace dsp
