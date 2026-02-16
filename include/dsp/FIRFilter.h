#pragma once

#include "Context.h"
#include "Delay.h"
#include "FastMath.h"
#include "Windows.h"

/*
    x0 |x1 x2 x3 x4| x5 x6 x7 x8 |x9
       |           |             |
 b6 b5 |b4 b3 b2 b1| b0 0  0  0  |
       |           |             |
 b7 b6 |b5 b4 b3 b2| b1 b0 0  0  |
       |           |             |
    b7 |b6 b5 b4 b3| b2 b1 b0 0  |
       |           |             |
       |b7 b6 b5 b4| b3 b2 b1 b0 |

*/

namespace dsp
{
inline namespace DSP_ARCH_NAMESPACE
{

template <typename T, size_t Order> class FIRFilter
{
  public:
    static constexpr auto kPad    = DSP_PADDING_VEC_SIZE / sizeof(T);
    static constexpr auto kNCoeff = Order + 1;

    template <typename T_ = T>
    class DL : public CopyDelayLine<T_, nextAlignedOffset(Order, kPad)>
    {
    };

    static constexpr auto kPaddedLength = kNCoeff + kPad * 2 - 1;

    FIRFilter() = default;
    FIRFilter(const std::array<T, kNCoeff> &b)
    {
        for (size_t k = 0; k < kNCoeff; ++k) {
            b_[kPaddedLength - kPad - k] = b[k];
        }
    }

    template <class Ctxt, class DL> void process(Ctxt ctxt, DL &delayline) const
    {
        auto x = ctxt.getInput();

        auto ctxtVec                = ctxt.vec();
        constexpr auto kIncrSize    = Ctxt::kIncrSize;
        constexpr auto kIncrSizeVec = decltype(ctxtVec)::kIncrSize;

        using simd           = decltype(ctxtVec.getInput());
        simd sums[kIncrSize] = {};

        for (size_t delay = 0; delay < kNCoeff + kIncrSizeVec - 1;
             delay += kIncrSizeVec) {
            simd x0;
            if (delay == 0) {
                if constexpr (Ctxt::kUseVec) {
                    x0 = ctxt.getInput();
                } else {
                    // trick to get input as vectorized (eg. {x,0,0,0})
                    T xtmp[decltype(ctxtVec)::kIncrSize]{};
                    ctxt.store(*xtmp, ctxt.getInput());
                    x0 = ctxtVec.load(*xtmp);
                }
            } else {
                x0 = delayline.read(ctxtVec, static_cast<int>(delay));
            }

            for (size_t n = 0; n < kIncrSize; ++n) {
                auto b = ctxtVec.load(b_[kPaddedLength - kPad - (delay + n)]);
                sums[n] += x0 * b;
            }
        }

        delayline.write(ctxt, x);

        T out[Ctxt::kIncrSize];
        for (size_t n = 0; n < Ctxt::kIncrSize; ++n) {
            out[n] = reduce<kTypeWidth<T>>(sums[n]);
        }
        ctxt.setOutput(ctxt.load(*out));
    }

  private:
    std::array<T, kPaddedLength> b_ = {};
};

template <typename T, size_t Order, size_t M,
          class Window = windows::Kaiser<140>>
class FIRDecimate
{
  public:
    static constexpr auto kPad    = DSP_PADDING_VEC_SIZE / sizeof(T);
    static constexpr auto kNCoeff = (Order + 1) * M;

    template <int Offset = 0>
    class DL : public DelayLine<nextAlignedOffset(kNCoeff - 1, kPad), Offset>
    {
    };

    static constexpr auto kPaddedLength = kNCoeff + kPad * 2 - 1;

    FIRDecimate(baseType<T> cutoff = 1)
    {
        auto freq = cutoff / M;
        auto mid  = (kNCoeff - 1) * baseType<T>(0.5);
        /* generate windowed sinc */
        for (size_t n = 0; n < kNCoeff; ++n) {
            auto fn = static_cast<baseType<T>>(n);
            b_[kPaddedLength - kPad - n] =
                Window::generate((fn - mid) / (mid)) * sinc((fn - mid) * freq) *
                freq;
        }

        /* scale */
        auto sum = load(T(0));
        for (auto &b : b_) {
            sum += load(b);
        }
        auto scale = T(1) / sum;
        for (auto &b : b_) {
            b = b * scale;
        }
    }

    const auto &getCoeffs() const { return b_; }

    template <class CtxtIn, class CtxtOut, class DL>
    int decimate(CtxtIn cin, CtxtOut &cout, DL &delayline, int decimateId) const
    {
        static_assert(!CtxtOut::kUseVec);
        auto decimatedBlockSize = 0;
        auto id                 = static_cast<size_t>(decimateId);
        auto coutCopy           = cout;

        CTXTRUN(cin)
        {
            constexpr auto kInIncrSize = decltype(cin)::kIncrSize;
            auto x                     = cin.getInput();

            size_t xOffset = (M - id) % M;
            while (xOffset < kInIncrSize) {
                decltype(x) sum = {};

                for (size_t delay = 0;
                     delay < kNCoeff + kInIncrSize - 1 - xOffset;
                     delay += kInIncrSize) {
                    auto x0 =
                        delay == 0
                            ? x
                            : delayline.read(cin, static_cast<int>(delay));
                    const auto b =
                        cin.load(b_[kPaddedLength - kPad - delay - xOffset]);

                    sum += x0 * b;
                }

                auto xout = dsp::reduce<kTypeWidth<T>>(sum);
                coutCopy.setOutput(xout);

                ++decimatedBlockSize;
                xOffset += M;
                coutCopy.next();
            }

            id = (id + kInIncrSize) % M;
            delayline.write(cin, x);
        };

        cout.setBlockSize(decimatedBlockSize);
        return static_cast<int>(id);
    }

  private:
    std::array<T, kPaddedLength> b_ = {};
};

/*
    x0 |x1 x2 x3 x4| x5 x6 x7 x8 |x9
    0  |0  0  b6 b3| b0 0  0  0  |
    0  |0  0  b7 b4| b1 0  0  0  |
    0  |0  0  b8 b5| b2 0  0  0  |
    0  |0  0  0  b6| b3 0  0  0  |
    0  |0  0  0  b7| b4 0  0  0  |
    0  |0  0  0  b8| b5 0  0  0  |
    0  |0  0  0  0 | b6 0  0  0  |
    0  |0  0  0  0 | b7 0  0  0  |


 */
template <typename T, size_t Order, size_t L,
          class Window = windows::Kaiser<140>>
class FIRInterpolate
{
  public:
    static constexpr auto kPad    = DSP_PADDING_VEC_SIZE / sizeof(T);
    static constexpr auto kNCoeff = (Order + 1);

    template <size_t Offset = 0>
    class DL : public DelayLine<nextAlignedOffset(kNCoeff, kPad), Offset>
    {
    };

    static constexpr auto kPaddedLength = kNCoeff + kPad * 2 - 1;

    FIRInterpolate(baseType<T> cutoff = 1)
    {
        auto freq = cutoff / L;
        auto mid  = (kNCoeff * L - 1) / 2.f;
        /* generate windowed sinc */
        for (size_t l = 0; l < L; ++l) {
            for (size_t n = 0; n < kNCoeff; ++n) {
                auto k  = l + n * L;
                auto fk = static_cast<baseType<T>>(k);
                b_[l][kPaddedLength - kPad - n] =
                    Window::generate((fk - mid) / (mid)) *
                    sinc((fk - mid) * freq) * cutoff;
            }
        }

        // scale
        auto sum = load(T(0));
        for (auto &b : b_)
            for (auto &bl : b) sum += bl;
        auto scale = T(L) / sum;
        for (auto &b : b_)
            for (auto &bl : b) bl = bl * scale;
    }

    template <class CtxtIn, class CtxtOut, class DL>
    int interpolate(CtxtIn cin, CtxtOut cout, DL &delayline,
                    int interpolateId) const
    {
        static_assert(CtxtIn::kIncrSize == 1);
        static_assert(CtxtOut::kIncrSize == 1);
        auto id = static_cast<size_t>(interpolateId);

        CTXTRUN(cout)
        {
            if (id == 0) {
                auto x = cin.getInput();
                delayline.writeSafe(cin, x);
                cin.next();
            }

            auto vcin                  = cin.vec();
            constexpr auto kInIncrSize = decltype(vcin)::kIncrSize;

            using outType = decltype(cout.getInput());
            using inType  = decltype(vcin.getInput());
            inType sum{};

            for (size_t delay = kInIncrSize; delay < kNCoeff + kInIncrSize;
                 delay += kInIncrSize) {
                auto x0 = delayline.read(vcin, static_cast<int>(delay));
                const auto b =
                    vcin.load(b_[id][kPaddedLength - kPad - (delay - 1)]);

                sum += x0 * b;
            }

            outType xout{};
            xout = dsp::reduce<kTypeWidth<T>>(sum);

            cout.setOutput(xout);

            id = (id + 1) % L;
        };
        return static_cast<int>(id);
    }

  private:
    T b_[L][kPaddedLength]{};
};

} // namespace DSP_ARCH_NAMESPACE
} // namespace dsp
