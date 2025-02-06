#pragma once

#include "Context.h"
#include "Delay.h"
#include "Utils.h"

namespace dsp
{

// https://www.dafx.de/paper-archive/2018/papers/DAFx2018_paper_9.pdf

template <int L> class TapePosition
{
  public:
    TapePosition() : tapePosBuf_{0, Unity} {}
    using position_t = int;

    static constexpr auto MinSize = L;
    static constexpr auto Size    = nextPow2(MinSize);
    static constexpr auto Mask    = Size - 1;

    static constexpr auto Q           = 30;
    static constexpr position_t Unity = (1 << Q);

    static constexpr position_t convertSpeed(float speed)
    {
        return Unity * speed;
    }

    void move(position_t speed)
    {
        tapePos_ += speed;

        tapePosBuf_[nwrite_] = tapePos_;
        nwrite_              = (nwrite_ + 1) & Mask;
    }

    int getNwrite() const { return nwrite_; }
    position_t get() const { return tapePos_; }

    position_t at(int index) const { return tapePosBuf_[index & Mask]; }

    static bool greaterThan(position_t x1, position_t x2)
    {
        unsigned int tmp = -x2;
        position_t diff  = x1 + tmp;
        return diff > 0;
    }

  private:
    int nwrite_{2};
    int count_{0};
    position_t tapePos_{Unity};
    position_t tapePosBuf_[Size];
};

template <class Tap = TapLin<1>> class TapTape
{
  public:
    template <class Ctxt, class DL, class TapePos>
    auto read(Ctxt c, const DL &delayline, const TapePos &tapePos)
    {
        // the tape position we need to read
        auto readPosition = tapePos.get() - TapePos::Unity;
        using position_t  = decltype(readPosition);

        // binary search of the position;
        auto nsearch = 2;
        size_t nread = nread_;
        for (position_t searchPosition;
             searchPosition = tapePos.at(nread + nsearch),
             TapePos::greaterThan(readPosition, searchPosition);
             nread += nsearch, nsearch *= 2) {

            // assert(nread + nsearch < tapePos.nwrite_);
        }

        while (nsearch > 1) {
            nsearch /= 2;
            auto searchPosition = tapePos.at(nread + nsearch);
            if (TapePos::greaterThan(readPosition, searchPosition)) {
                nread += nsearch;
            }
        }

        nread_ = nread & TapePos::Mask;

        auto foundPos0 = tapePos.at(nread_);
        auto foundPos1 = tapePos.at(nread_ + 1);

        // determine delay;
        int delay   = (tapePos.getNwrite() - nread_ + 1) & TapePos::Mask;
        auto fdelay = static_cast<float>(foundPos1 - readPosition) /
                      (foundPos1 - foundPos0);

        Tap tap;
        tap.setDelay({delay}, {fdelay});
        return tap.read(c, delayline);
    }

  private:
    // index of last read position
    int nread_{0};
};
} // namespace dsp
