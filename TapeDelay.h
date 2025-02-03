#pragma once

#include "Delay.h"
#include "Utils.h"

namespace dsp
{

// https://www.dafx.de/paper-archive/2018/papers/DAFx2018_paper_9.pdf

class TapTape;

template <int L> class TapePosition
{
  public:
    friend TapTape;

    TapePosition() : tapePosBuf_{0, Unity} {}
    using position_t = int;

    static constexpr auto MinSize = L;
    static constexpr auto Size    = nextPow2(MinSize);
    static constexpr auto Mask    = Size - 1;

    static constexpr auto Q           = 30;
    static constexpr position_t Unity = (1 << Q);

    void move(float speed)
    {
        position_t ispeed = Unity * speed;
        tapePos_ += ispeed;

        tapePosBuf_[nwrite_] = tapePos_;
        ++nwrite_;
    }

    position_t get() const { return tapePos_; }

    position_t at(std::size_t index) const { return tapePosBuf_[index & Mask]; }

    static bool greaterThan(position_t x1, position_t x2)
    {
        auto diff = x1 - x2;
        return diff > 0;
    }

  private:
    std::size_t nwrite_{2};
    int count_{0};
    position_t tapePos_{Unity};
    position_t tapePosBuf_[Size];
};

class TapTape
{
  public:
    template <class Ctxt, class DL, class TapePos>
    auto read(Ctxt c, DL &delayline, TapePos tapePos)
    {
        // the tape position we need to read
        auto readPosition = tapePos.get() - TapePos::Unity;
        using position_t  = decltype(readPosition);

        // binary search of the position;
        auto nsearch = 1;
        for (position_t searchPosition;
             searchPosition = tapePos.at(nread_ + nsearch),
             TapePos::greaterThan(readPosition, searchPosition);
             nread_ += nsearch, nsearch *= 2) {

            assert(nread_ + nsearch < tapePos.nwrite_);
        }

        while (nsearch > 1) {
            nsearch /= 2;
            auto searchPosition = tapePos.at(nread_ + nsearch);
            if (TapePos::greaterThan(readPosition, searchPosition)) {
                nread_ += nsearch;
            }
        }

        nread_ &= TapePos::Mask;

        auto foundPos0 = tapePos.at(nread_);
        auto foundPos1 = tapePos.at(nread_ + 1);

        // determine delay;
        auto delay  = (tapePos.nwrite_ - nread_);
        auto fdelay = static_cast<float>(readPosition - foundPos0) /
                      (foundPos1 - foundPos0);

        auto x0 = delayline.read(c, delay);
        auto x1 = delayline.read(c, delay + 1);
        for (int k = 0; k < x0.size(); ++k) {
            for (int i = 0; i < x0[0].size(); ++i) {
                x1[k][i] += fdelay * (x0[k][i] - x1[k][i]);
            }
        }
        return x1;
    }

  private:
    // index of last read position
    std::size_t nread_{0};
};
} // namespace dsp
