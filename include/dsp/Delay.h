#pragma once

#include "Context.h"
#include "FastMath.h"
#include "Utils.h"
#include <cassert>

namespace dsp
{

class BufferOffset
{
  public:
    int add(int length)
    {
        const auto offset = value_;
        value_ += length + 1;
        return offset;
    }

  private:
    int value_{0};
};

static BufferOffset bufferOffset;

template <size_t L = 0, int Off = 0> class DelayLine;

template <> class DelayLine<>
{
  public:
    /*
     * Implement a delay line
     * templates:
     *      N: vector size
     *      Length: maximum length of the delay line
     */

    static constexpr size_t kLength     = 10000000;
    static constexpr size_t kNextOffset = 0;

    template <int> using WithOffset = DelayLine;

    DelayLine(int length) : offset_(bufferOffset.add(length)) {}

    /* write value inside delay line */
    template <class Ctxt, class In> void write(Ctxt c, const In &x) const
    {
        c.write(offset_, x);
    }
    template <class Ctxt, class In> void writeSafe(Ctxt c, const In &x) const
    {
        c.writeSafe(offset_, x);
    }

    /* read value at delay id */
    template <class Ctxt> const auto &read(Ctxt c, int id) const
    {
        return c.read(offset_ + id);
    }

  private:
    const int offset_;
};

/* delayline with compile time length */
template <size_t L, int Off> class DelayLine
{
    /* find the next vector aligned offset */
  public:
    static constexpr auto kLength = L;
    static constexpr auto kOffset = Off;
    static constexpr auto kNextOffset =
        kOffset + nextAlignedOffset(kLength, DSP_MAX_VEC_SIZE / 4);

    template <int O> using WithOffset = DelayLine<L, O>;

    constexpr DelayLine() = default;

    /* write value inside delay line */
    template <class Ctxt, class In> void write(Ctxt c, const In &x) const
    {
        c.write(kOffset, x);
    }
    template <class Ctxt, class In> void writeSafe(Ctxt c, const In &x) const
    {
        c.writeSafe(kOffset, x);
    }

    /* read value at delay id */
    template <class Ctxt> [[nodiscard]] auto read(Ctxt c, int id) const
    {
        return c.read(kOffset + id);
    }

    /* read value at tail */
    template <class Ctxt> auto tail(Ctxt c) const { return read(c, kLength); }
};

template <typename T, size_t L = 1, int Off = 0> class CopyDelayLine
{
  public:
    static_assert(L > 0, "Length must be bigger than 0");
    static constexpr auto kLength     = L;
    static constexpr auto kOffset     = Off;
    static constexpr auto kNextOffset = kOffset;

    template <int O> using WithOffset = CopyDelayLine<T, L, O>;

    template <class Ctxt, typename V> void write(Ctxt c, V x)
    {
        (void)c;

        constexpr auto kVecSize = Ctxt::kIncrSize;
        auto lastNonVector      = kLength % kVecSize;

        /* first shift non vector aligned float */
        for (size_t j = 0; j < lastNonVector; ++j) mem_[j] = mem_[j + kVecSize];

        /* shift rest */
        for (size_t j = lastNonVector; j < kLength - kVecSize; j += kVecSize)
            c.store(mem_[j], c.load(mem_[j + kVecSize]));

        /* copy last value */
        c.store(mem_[kLength - kVecSize], x);
    }
    template <class Ctxt, typename V> void writeSafe(Ctxt c, V x)
    {
        write(c, x);
    }

    template <class Ctxt> [[nodiscard]] auto read(Ctxt c, int i) const
    {
        auto si = static_cast<size_t>(i);

        assert(0 <= i && si <= kLength);
        return c.load(mem_[kLength - si]);
    }

    template <class Ctxt> [[nodiscard]] auto tail(Ctxt c) const
    {
        return read(c, kLength);
    }

  private:
    T mem_[kLength]{};
};

template <class DL, class DLi, int Off = 0>
class NestedDelayLine : public DL::template WithOffset<Off>
{
    using Outer = typename DL::template WithOffset<Off>;
    using Inner = typename DLi::template WithOffset<Outer::NextOffset>;

  public:
    static constexpr auto kLength     = Outer::Length + Inner::Length;
    static constexpr auto kOffset     = Off;
    static constexpr auto kNextOffset = Inner::NextOffset;

    template <int O> using WithOffset = NestedDelayLine<DL, DLi, O>;

    NestedDelayLine() = default;
    NestedDelayLine(Outer &&dl, Inner &&dli) : Outer(dl), inner_(dli) {}

    Inner &getInner() { return inner_; }

  private:
    Inner inner_;
};

template <class DL, size_t Nm, int Off = 0>
class ArrayDelayLine
    : protected std::array<typename DL::template WithOffset<Off>, Nm>
{
    using Base = typename DL::template WithOffset<Off>;

  public:
    static constexpr auto kLength     = DL::Length * Nm;
    static constexpr auto kOffset     = Off;
    static constexpr auto kOneOffset  = (Base::NextOffset - kOffset);
    static constexpr auto kNextOffset = kOneOffset * Nm;

    template <int O> using WithOffset = ArrayDelayLine<DL, Nm, O>;

    template <class Ctxt> auto &get(Ctxt &c, int i)
    {
        c.nextBufId(kOneOffset * i);
        return (*this)[i];
    }
};

#define nextTo(delayline) decltype(delayline)::kNextOffset

/*
 * Taps: different classes that can read into a delay line with different
 * position and interpolation In each class DelayLineType<L> define the type of
 * the base DelayLine that can be read by the tap
 */

struct TapTail {
    /* Tap that reads at the end of a delayline */

    template <class Ctxt, class DL>
    [[nodiscard]] auto read(Ctxt c, const DL &delayline) const
    {
        return delayline.tail(c);
    }
};

/* help function */
template <class Ctxt, class DL, class T>
void basefixread(Ctxt /*c*/, const DL & /*delayline*/, T & /*x*/, int /*i*/)
{
}
template <class Ctxt, class DL, class T, int D, int... Ds>
void basefixread(Ctxt c, const DL &delayline, T &x, int i)
{
    static_assert(D <= DL::kLength - Ctxt::Type::kWidth + 1,
                  "tap delay length is bigger than delay line");
    auto val = delayline.read(c, D);
    x[i]     = val[i];
    basefixread<Ctxt, DL, T, Ds...>(c, delayline, x, ++i);
}
template <int D = 1, int... Ds> struct TapFix : public TapFix<Ds...> {
    template <class Ctxt, class DL> auto read(Ctxt c, const DL &delayline) const
    {
        typename Ctxt::Type x = {0};
        basefixread<Ctxt, DL, decltype(x), D, Ds...>(c, delayline, x, 0);
        return c.load(x);
    }
};
template <int D> struct TapFix<D> {
    /* Tap that reads at a fix point (delay D) in the delay line */

    template <class Ctxt, class DL> auto read(Ctxt c, const DL &delayline) const
    {
        static_assert(D > 0 && D <= DL::kLength - Ctxt::kIncrSize + 1,
                      "tap delay length is bigger than delay line");
        return delayline.read(c, D);
    }
};
template <typename T> class TapNoInterp
{
    /* Tap that can read at any point in de delay line, without interpolation */
  public:
    using IntT = intType<T>;

    TapNoInterp() = default;
    TapNoInterp(IntT id) : id_(id) {}

    void setDelay(int i, int id) { id_[i] = id; }
    void setDelay(const IntT &id) { id_ = id; }

    template <class Ctxt, class DL>
    [[nodiscard]] auto read(Ctxt c, const DL &delayline) const
    {
        constexpr auto kWidth      = kTypeWidth<typename Ctxt::Type>;
        constexpr auto kDelayWidth = kTypeWidth<IntT>;

        constexpr int kMinDelay = Ctxt::kIncrSize;
        (void)kMinDelay;
        assert(all(load(id_) >= kMinDelay &&
                   load(id_) <= static_cast<int>(DL::kLength - kMinDelay + 1)));

        if constexpr (kDelayWidth > 1) {
            typename Ctxt::Type x;
            for (size_t i = 0; i < kWidth; ++i) {
                x[i] = delayline.read(c, id_[i % kDelayWidth])[i];
            }
            return c.load(x);
        } else {
            return delayline.read(c, id_);
        }
    }

  protected:
    IntT id_;
};

template <typename T> class TapLin : public TapNoInterp<T>
{
    /* Tap that can read at any point in de delay line, with linear
     * interpolation */
  public:
    using IntT = intType<T>;

    TapLin() = default;
    TapLin(T d) { setDelay(d); }

    void setDelay(const T &d)
    {
        IntT id;
        store(id, toInt(load(d)));
        store(fd_, load(d) - id);
        TapNoInterp<T>::setDelay(id);
    }

    void setDelay(const IntT &d, const T &fd)
    {
        TapNoInterp<T>::setDelay(d);

        auto vfd = load(fd);
        (void)vfd;
        assert(all(vfd <= T(1) && vfd >= T(0)));

        fd_ = fd;
    }

    template <class Ctxt, class DL> auto read(Ctxt c, DL &delayline) const
    {
        auto id = load(TapNoInterp<T>::id_);

        assert(all(id > 0 &&
                   id <= static_cast<int>(DL::kLength - Ctxt::kIncrSize)));

        auto x0 = TapNoInterp<T>::read(c, delayline);
        auto x1 = TapNoInterp<T>{id + 1}.read(c, delayline);

        x0 += fd_ * (x1 - x0);
        return x0;
    }

  protected:
    T fd_;
};

template <typename T> class TapCubic : public TapLin<T>
{
    /* Tap that can read at any point in the delay line, with cubic
     * interpolation */

  public:
    TapCubic() = default;
    TapCubic(const T &delay) : TapLin<T>(delay) {}

    template <class Ctxt, class DL> auto read(Ctxt c, DL &delayline) const
    {
        auto id = load(TapNoInterp<T>::id_);
        auto fd = load(TapLin<T>::fd_);

        assert(any(id > 1 &&
                   id <= static_cast<int>(DL::kLength - Ctxt::kIncrSize - 2)));

        auto xm1 = TapNoInterp<T>{id - 1}.read(c, delayline);
        auto x0  = TapNoInterp<T>::read(c, delayline);
        auto x1  = TapNoInterp<T>{id + 1}.read(c, delayline);
        auto x2  = TapNoInterp<T>{id + 2}.read(c, delayline);

        return hermite(xm1, x0, x1, x2, fd);
    }
};
} // namespace dsp
