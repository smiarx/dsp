#pragma once

#include "simd/multi.h"
#include <limits>
#include <random>

namespace dsp
{

template <typename T> class INoise
{
    using iT                          = MultiVal<int, kTypeWidth<T>>;
    static constexpr int kGenerator[] = {1296462733, 1487623987, 848278349,
                                         1987647829, 1837654627, 963782763,
                                         1492758273, 1746273223};

  public:
    INoise() : INoise(std::random_device{}()) {}
    INoise(int seed)
    {
        auto rnd = seed;
        for (size_t i = 0; i < kTypeWidth<iT>; ++i) {
            seed_[i] = rnd;
            rnd      = (rnd + seed) * 1103515245;
        }
    }

    auto process()
    {
        state_ = seed_ + state_ * iT::load(kGenerator);
        return load(state_);
    }

  private:
    iT seed_;
    iT state_{0};
};

template <typename T> class Noise : public INoise<T>
{
  public:
    Noise(int seed) : INoise<T>(seed) {}
    Noise() = default;
    auto process()
    {
        using type  = baseType<T>;
        using itype = intType<T>;

        auto iy = INoise<T>::process();
        auto y  = toFloat<type>(iy) /
                 static_cast<type>(std::numeric_limits<baseType<itype>>::max());
        return y;
    }
};

} // namespace dsp
