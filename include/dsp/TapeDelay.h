#pragma once

#include "Delay.h"
#include "Utils.h"

namespace dsp
{

// https://www.dafx.de/paper-archive/2018/papers/DAFx2018_paper_9.pdf

template <size_t L> class TapePosition
{
  public:
    TapePosition() : tapePosBuf_{0, Unity} {}
    using position_t = int;

    static constexpr auto MinSize = L;
    static constexpr auto Size    = nextPow2(MinSize);
    static constexpr auto Mask    = Size - 1;

    static constexpr auto Q           = 29;
    static constexpr position_t Unity = (1 << Q);

    static constexpr position_t convertSpeed(float speed)
    {
        return static_cast<position_t>(Unity * speed);
    }

    void move(position_t speed)
    {
        tapePos_ += speed;

        tapePosBuf_[nwrite_] = tapePos_;
        nwrite_              = (nwrite_ + 1) & Mask;
    }

    [[nodiscard]] size_t getNwrite() const { return nwrite_; }
    [[nodiscard]] position_t get() const { return tapePos_; }

    [[nodiscard]] position_t at(size_t index) const
    {
        return tapePosBuf_[index & Mask];
    }

    static bool greaterThan(position_t x1, position_t x2)
    {
        auto tmp = static_cast<unsigned int>(-x2);
        auto diff =
            static_cast<position_t>(static_cast<unsigned int>(x1) + tmp);
        return diff > 0;
    }

  private:
    size_t nwrite_{2};
    int count_{0};
    position_t tapePos_{Unity};
    position_t tapePosBuf_[Size];
};

template <class Tap = TapLin<1>> class TapTape
{
  public:
    enum Direction {
        Normal,
        Reverse,
    };

    // search given position in tape
    template <class TapePos>
    auto search(const TapePos &tapePos,
                typename TapePos::position_t distance = TapePos::Unity)
    {
        // the tape position we need to read
        auto readPosition = tapePos.get() - distance;
        using position_t  = decltype(readPosition);

        // binary search of the position;
        constexpr Direction D = Reverse;
        size_t nsearch        = 2;
        auto nread            = static_cast<size_t>(tapePos.getNwrite());
        for (position_t searchPosition;
             searchPosition = tapePos.at(add<D>(nread, nsearch)),
             greaterThan<D, TapePos>(readPosition, searchPosition);
             nread = add<D>(nread, nsearch), nsearch *= 2) {
        }

        while (nsearch > 1) {
            nsearch /= 2;
            auto nreadnew       = add<D>(nread, nsearch);
            auto searchPosition = tapePos.at(nreadnew);
            if (greaterThan<D, TapePos>(readPosition, searchPosition)) {
                nread = nreadnew;
            }
        }

        nread &= TapePos::Mask;
        nread_ = nread - 1;
    }

    template <Direction D = Normal, class Ctxt, class DL, class TapePos>
    auto read(Ctxt c, const DL &delayline, const TapePos &tapePos,
              typename TapePos::position_t distance = TapePos::Unity)
    {
        // the tape position we need to read
        auto readPosition = tapePos.get() - distance;
        using position_t  = decltype(readPosition);

        // binary search of the position;
        size_t nsearch = 2;
        size_t nread   = nread_;
        for (position_t searchPosition;
             searchPosition = tapePos.at(add<D>(nread, nsearch)),
             greaterThan<D, TapePos>(readPosition, searchPosition);
             nread = add<D>(nread, nsearch), nsearch *= 2) {
        }

        while (nsearch > 1) {
            nsearch /= 2;
            auto nreadnew       = add<D>(nread, nsearch);
            auto searchPosition = tapePos.at(nreadnew);
            if (greaterThan<D, TapePos>(readPosition, searchPosition)) {
                nread = nreadnew;
            }
        }

        nread &= TapePos::Mask;
        nread_ = nread;

        if constexpr (D == Reverse) {
            nread = (nread - 1) & TapePos::Mask;
        }

        auto foundPos0 = tapePos.at(nread);
        auto foundPos1 = tapePos.at((nread + 1) & TapePos::Mask);

        // determine delay;
        int delay   = (tapePos.getNwrite() - nread_ + 1) & TapePos::Mask;
        auto fdelay = static_cast<float>(foundPos1 - readPosition) /
                      static_cast<float>(foundPos1 - foundPos0);

        Tap tap;
        tap.setDelay({delay}, {fdelay});
        //(void)fdelay;
        // tap.setDelay({delay,delay}, {fdelay,fdelay});
        return tap.read(c, delayline);
    }

    template <class TapePos> void reset(TapePos tapePos)
    {
        nread_ = tapePos.getNwrite();
    }

  private:
    template <Direction D> size_t add(size_t n1, size_t n2)
    {
        if constexpr (D == Normal) {
            return n1 + n2;
        } else {
            return n1 - n2;
        }
    }

    template <Direction D, class TapePos>
    bool greaterThan(typename TapePos::position_t p1,
                     typename TapePos::position_t p2)
    {
        if constexpr (D == Normal) {
            return TapePos::greaterThan(p1, p2);
        } else {
            return TapePos::greaterThan(p2, p1);
        }
    }

    // index of last read position
    size_t nread_{0};
};

} // namespace dsp
