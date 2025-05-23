#pragma once

#include "Context.h"
#include "FastMath.h"
#include "Signal.h"
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
        c.getBuffer().write(offset_, x);
    }

    /* read value at delay id */
    template <class Ctxt> const auto &read(Ctxt c, int id) const
    {
        return c.getBuffer().read(offset_ + id);
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
        kOffset + nextAlignedOffset(kLength, SIMDSIZE);

    template <int O> using WithOffset = DelayLine<L, O>;

    constexpr DelayLine() = default;

    /* write value inside delay line */
    template <class Ctxt, class In> void write(Ctxt c, const In &x) const
    {
        c.getBuffer().write(kOffset, x);
    }

    /* read value at delay id */
    template <class Ctxt> [[nodiscard]] const auto &read(Ctxt c, int id) const
    {
        return c.read(kOffset + id);
    }

    /* read value at tail */
    template <class Ctxt> const auto &tail(Ctxt c) const
    {
        return read(c, kLength);
    }
};

template <size_t N, size_t L = 1, int Off = 0> class CopyDelayLine
{
  public:
    static_assert(L > 0, "Length must be bigger than 0");
    static constexpr auto kLength     = L;
    static constexpr auto kOffset     = Off;
    static constexpr auto kNextOffset = kOffset;

    template <int O> using WithOffset = CopyDelayLine<N, L, O>;

    template <class Ctxt> void write(Ctxt c, const typename Ctxt::Type &x)
    {
        (void)c;

        constexpr auto kVecSize = Ctxt::kVecSize;
        constexpr auto kIsVec   = Ctxt::kIsUsingVector;
        auto lastNonVector      = kLength % kVecSize;

        /* first shift non vector aligned float */
        for (size_t j = 0; j < lastNonVector; ++j) mem_[j] = mem_[j + kVecSize];

        /* shift rest */
        for (size_t j = lastNonVector; j < kLength - kVecSize; j += kVecSize)
            mem_[j].template toSignal<kIsVec>() =
                mem_[j + kVecSize].template toSignal<kIsVec>();

        /* copy last value */
        mem_[kLength - kVecSize].template toSignal<kIsVec>() = x;
    }

    template <class Ctxt> [[nodiscard]] const auto &read(Ctxt c, int i) const
    {
        (void)c;
        auto si = static_cast<size_t>(i);

        assert(si <= kLength);
        const auto &x = mem_[kLength - si];
        if constexpr (Ctxt::kIsUsingVector) {
            return x.toVector();
        } else {
            return x.toScalar();
        }
    }

    template <class Ctxt> [[nodiscard]] const auto &tail(Ctxt c) const
    {
        return read(c, kLength);
    }

  private:
    fSample<N> mem_[kLength]{};
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
    [[nodiscard]] const auto &read(Ctxt c, const DL &delayline) const
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
    static_assert(D <= DL::Length - Ctxt::VecSize + 1,
                  "tap delay length is bigger than delay line");
    auto &val = delayline.read(c, D);
    for (size_t k = 0; k < Ctxt::VecSize; ++k) x[k][i] = val[k][i];
    basefixread<Ctxt, DL, T, Ds...>(c, delayline, x, ++i);
}
template <int D = 1, int... Ds> struct TapFix : public TapFix<Ds...> {
    template <class Ctxt, class DL> auto read(Ctxt c, const DL &delayline) const
    {
        typename Ctxt::Type x = {0};
        basefixread<Ctxt, DL, decltype(x), D, Ds...>(c, delayline, x, 0);
        return x;
    }
};
template <int D> struct TapFix<D> {
    /* Tap that reads at a fix point (delay D) in the delay line */

    template <class Ctxt, class DL>
    const auto &read(Ctxt c, const DL &delayline) const
    {
        static_assert(D <= DL::Length - Ctxt::VecSize + 1,
                      "tap delay length is bigger than delay line");
        return delayline.read(c, D);
    }
};
template <size_t N> class TapNoInterp
{
    /* Tap that can read at any point in de delay line, without interpolation */
  public:
    void setDelay(int i, int id) { id_[i] = id; }
    void setDelay(iData<N> id) { id_ = id; }

    template <class Ctxt, class DL>
    [[nodiscard]] auto read(Ctxt c, const DL &delayline) const
    {
        typename Ctxt::Type x;
        for (size_t i = 0; i < N; i++) {
            assert(id_[i] <=
                   static_cast<int>(DL::kLength - Ctxt::kVecSize + 1));
            auto &val = delayline.read(c, id_[i]);
            for (size_t k = 0; k < Ctxt::kVecSize; ++k) x[k][i] = val[k][i];
        }
        return x;
    }

  protected:
    iData<N> id_;
};

template <size_t N> class TapLin : public TapNoInterp<N>
{
    /* Tap that can read at any point in de delay line, with linear
     * interpolation */
  public:
    void setDelay(fData<N> d)
    {
        iData<N> id;
        for (size_t i = 0; i < N; ++i) {
            id[i]  = static_cast<int>(d[i]);
            fd_[i] = d[i] - static_cast<float>(id[i]);
        }
        TapNoInterp<N>::setDelay(id);
    }

    void setDelay(iData<N> d, fData<N> fd)
    {
        TapNoInterp<N>::setDelay(d);

        arrayFor(fd, i) { assert(fd[i] <= 1.f && fd[i] >= 0.f); }
        fd_ = fd;
    }

    template <class Ctxt, class DL> auto read(Ctxt c, DL &delayline) const
    {
        typename Ctxt::Type x0;
        typename Ctxt::Type x1;
        auto id = TapNoInterp<N>::id_;

        if constexpr (N == 1) {
            assert(id[0] <= static_cast<int>(DL::Length - Ctxt::VecSize - 1));
            x0 = delayline.read(c, id[0]);
            x1 = delayline.read(c, id[0] + 1);
        } else {
            arrayFor(x0[0], i)
            {
                assert(id[i % N] <= DL::Length - Ctxt::VecSize - 1);
                auto &val0 = delayline.read(c, id[i % N]);
                for (size_t k = 0; k < Ctxt::VecSize; ++k)
                    x0[k][i] = val0[k][i];
                auto &val1 = delayline.read(c, id[i % N] + 1);
                for (size_t k = 0; k < Ctxt::VecSize; ++k)
                    x1[k][i] = val1[k][i];
            }
        }

        inFor(x0, k, i) { x0[k][i] += fd_[i % N] * (x1[k][i] - x0[k][i]); }
        return x0;
    }

  protected:
    fData<N> fd_;
};

template <size_t N> class TapCubic : public TapLin<N>
{
    /* Tap that can read at any point in the delay line, with cubic
     * interpolation */

  public:
    template <class Ctxt, class DL> auto read(Ctxt c, DL &delayline) const
    {
        typename Ctxt::Type y;
        typename Ctxt::Type xm1;
        typename Ctxt::Type x0;
        typename Ctxt::Type x1;
        typename Ctxt::Type x2;
        auto id = TapNoInterp<N>::id_;
        auto fd = TapLin<N>::fd_;

        if constexpr (N == 1) {
            assert(id[0] <= DL::Length - Ctxt::VecSize - 2);
            xm1 = delayline.read(c, id[0] - 1);
            x0  = delayline.read(c, id[0]);
            x1  = delayline.read(c, id[0] + 1);
            x2  = delayline.read(c, id[0] + 2);
        } else {
            arrayFor(x0[0], i)
            {
                assert(id[i % N] <=
                       static_cast<int>(DL::kLength - 2 - Ctxt::kVecSize));
                auto &valm1 = delayline.read(c, id[i % N] - 1);
                for (size_t k = 0; k < Ctxt::kVecSize; ++k)
                    xm1[k][i] = valm1[k][i];
                auto &val0 = delayline.read(c, id[i % N]);
                for (size_t k = 0; k < Ctxt::kVecSize; ++k)
                    x0[k][i] = val0[k][i];
                auto &val1 = delayline.read(c, id[i % N] + 1);
                for (size_t k = 0; k < Ctxt::kVecSize; ++k)
                    x1[k][i] = val1[k][i];
                auto &val2 = delayline.read(c, id[i % N] + 2);
                for (size_t k = 0; k < Ctxt::kVecSize; ++k)
                    x2[k][i] = val2[k][i];
            }
        }

        inFor(y, k, i)
        {
            y[k][i] =
                hermite(xm1[k][i], x0[k][i], x1[k][i], x2[k][i], fd[i % N]);
        }
        return y;
    }
};
} // namespace dsp
