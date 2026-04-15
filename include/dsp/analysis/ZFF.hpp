#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <dsp/cpu/infos.h>

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T, size_t MaxMarkers = 4> class ZFF
{
    // zero frequency filtering epoch analysis
  public:
    struct Marker {
        T position_;
        T period_;
    };

    ZFF() = default;
    ZFF(int windowSize) { setWindowSize(windowSize); }

    void setWindowSize(int windowSize)
    {
        windowSize_    = windowSize;
        invWindowSize_ = T(1) / static_cast<T>(windowSize);
    }

    void setLeak(T leak) { leak_ = leak; }

    constexpr auto getMaxMarkers() const { return MaxMarkers; }
    const auto &getMarkers() const { return markers_; }

    static constexpr bool defaultCallback(const Marker &) { return true; }

    template <class Ctxt>
    auto &getClosestMarker(Ctxt ctxt, T playRate = T(1))
    {
        auto foundId = markerId_;
        while (ctxt.idDiff(markers_[foundId].position_) -
                   markers_[foundId].period_ * (playRate - T(1)) <
               T(0)) {
            foundId = (foundId - 1 + MaxMarkers) % MaxMarkers;

            // going full circle
            if(foundId == markerId_) assert(false);
        }
        return markers_[foundId];
    }

    void addMarker(const Marker &newMarker)
    {
        // assert(newMarker.period_ > T(0));
        markerId_ = (markerId_ + 1) % MaxMarkers;
        nMarkers_ = std::min(nMarkers_ + 1, MaxMarkers);
        assert(newMarker.period_ >= T(0));
        markers_[markerId_] = newMarker;
    }

    const auto &getMarker(size_t id = 0) const
    {
        return markers_[(markerId_ + MaxMarkers - id) % MaxMarkers];
    }
    auto &getMarker(size_t id = 0) 
    {
        return markers_[(markerId_ + MaxMarkers - id) % MaxMarkers];
    }

    auto getNumMarkers() const { return nMarkers_; }

    void resetMarkers() { nMarkers_ = 0; }

    template <class Ctxt, class DL, class Callback = decltype(defaultCallback)>
    void process(Ctxt ctxt, DL delayline, Callback callback = defaultCallback)
    {
        auto x = ctxt.getInput();

        // differentiate
        auto xd = x - xPrev_;
        xPrev_  = x;

        // 2 stage leaky integrator
        y1_ += xd - leak_ * y1_;
        y2_ += y1_ - leak_ * y2_;

        // trend removal
        delayline.write(ctxt, y2_);
        auto oldY2 = delayline.read(ctxt, windowSize_);
        runningSum_ += y2_ - oldY2;

        auto trend = runningSum_ * invWindowSize_;
        auto zff   = y2_ - trend;

        // compute slope
        slope_ += T(0.03) * (zff - zffPrev_ - slope_);

        // zero crossing
        if (zffPrev_ < T(0) && zff >= T(0)) { // && slope_ > slopeThresh_) {
            T position;
            T denom = zff - zffPrev_;
            if (denom != 0) {
                // linear interpolation
                auto alpha = -zffPrev_ / denom;
                position   = static_cast<T>(ctxt.idDiff(1)) + alpha;
            } else
                position = static_cast<T>(ctxt.getId());

            Marker newMarker{position,
                             nMarkers_ > 0
                                 ? ctxt.modDiff(position, getMarker().position_)
                                 : T(0)};

            // we won't write marker before last marker
            if (newMarker.period_ >= T(0)) {
                if (callback(newMarker)) {
                    addMarker(newMarker);
                }
            }
        }

        zffPrev_ = zff;
    }

  private:
    T leak_{T(3e-3)};
    T slopeThresh_{T(0)};
    T slope_{};

    int windowSize_{};
    T invWindowSize_{};

    T xPrev_{};
    T y1_{}, y2_{};
    T zffPrev_{};

    T runningSum_{};

    size_t markerId_{};
    size_t nMarkers_{};
    std::array<Marker, MaxMarkers> markers_{};
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
