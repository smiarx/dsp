#include <catch2/catch_all.hpp>
#include <dsp/core/Buffer.hpp>
#include <dsp/math/FastMath.hpp>
#include <dsp/pitch/NSDF.hpp>

template <class N> static auto findPeakPeriod(N &nsdf)
{
    auto nPeaks       = nsdf.updatePeaks();
    const auto &peaks = nsdf.getPeaks();
    dsp::interp::PolyPoint<float> maxPeak{};
    for (size_t i = 0; i < nPeaks; ++i) {
        const auto &peak = peaks[i];
        if (peak.y > maxPeak.y) {
            maxPeak = peak;
        }
    }

    for (size_t i = 0; i < nPeaks; ++i) {
        const auto &peak = peaks[i];
        if (peak.y > maxPeak.y * 0.9f) {
            return peak.x;
        }
    }
    return maxPeak.x;
}

TEST_CASE("NSDF")
{
    dsp::NSDF<float, 100> nsdf;
    decltype(nsdf)::DL<0> dl;
    dsp::Buffer<float, nextTo(dl)> buffer;
    float bufferData[decltype(buffer)::kSize]{};
    buffer.setData(bufferData);

    SECTION("Autocorrelation")
    {
        float randx[] = {
            0.32165611783132775f, 0.682729830435953f,   0.46150737884763116f,
            0.5214082242322124f,  0.7437585981332573f,  0.47469329717700026f,
            0.7073313036108773f,  0.38635759075613973f, 0.977798859207111f,
            0.09479469049536682f, 0.6195859525937197f,  0.6004364269285843f,
            0.15858144060968682f, 0.8423231788235385f,  0.8790811185776595f,
            0.5670729010325121f,  0.29034657360327454f, 0.1770546143487549f,
            0.9298937730791831f,  0.46356077167956145f, 0.8857869724873033f,
            0.5521366506002903f,  0.6029347249386067f,  0.8847620676955268f,
            0.03214374284001398f, 0.3305275110901924f,  0.41704126780617623f,
            0.4554543413797768f,  0.18413704247453522f, 0.6918312611453639f,
            0.2609961823667565f,  0.8840458914352928f,  0.4958553317520731f,
            0.8578014309297725f,  0.3403940779238238f,  0.9651153921963733f,
            0.3790416337232465f,  0.32562060050080477f, 0.7376684560345416f,
            0.6140574184957287f,  0.7224684474742655f,  0.834637797931174f,
            0.2805682340539434f,  0.9167385552906636f,  0.3010393882345608f,
            0.03974294534681311f, 0.14692931368054585f, 0.578570259066359f,
            0.511275811827078f,   0.0800644325355293f,  0.24036853894700239f,
            0.7923200001279825f,  0.8611135200244278f,  0.3808157132341914f,
            0.42804526707331636f, 0.6063605913790103f,  0.5193444398959609f,
            0.5085452995726278f,  0.7213969048996709f,  0.6958171551643914f};
        constexpr auto kN = sizeof(randx) / sizeof(randx[0]);

        dsp::BufferContext ctxt(randx, kN, buffer);
        CTXTRUN(ctxt) { nsdf.process(ctxt, dl); };
    }

    constexpr auto kN = 300;
    float x[kN]{};

    SECTION("Integer Periods")
    {
        int period = GENERATE(5, 18, 30, 60);

        for (size_t n = 0; n < kN; ++n) {
            x[n] = dsp::sin(2 * dsp::constants<float>::pi /
                            static_cast<float>(period) * static_cast<float>(n));
        }

        dsp::BufferContext ctxt(x, kN, buffer);
        CTXTRUN(ctxt) { nsdf.process(ctxt, dl); };

        // multiples of period have nsdf of 1
        for (auto t = static_cast<size_t>(period);
             t < decltype(nsdf)::kTauLength; t += static_cast<size_t>(period)) {
            REQUIRE_THAT(nsdf.get(t), Catch::Matchers::WithinAbs(1.f, 1e-5f));
        }
        auto peakPeriod = findPeakPeriod(nsdf);
        REQUIRE_THAT(peakPeriod, Catch::Matchers::WithinAbs(period, 5e-2));
    }
    SECTION("Harmonics")
    {
        auto period = GENERATE(15.7872f, 30.3723f, 70.8827f);
        for (size_t n = 0; n < kN; ++n) {
            x[n] = dsp::sin(2 * dsp::constants<float>::pi / period *
                            static_cast<float>(n)) +
                   1.045f * dsp::sin(2 * 2 * dsp::constants<float>::pi /
                                     period * static_cast<float>(n)) +
                   0.875f * dsp::sin(3 * 2 * dsp::constants<float>::pi /
                                     period * static_cast<float>(n)) +
                   0.677f * dsp::sin(4 * 2 * dsp::constants<float>::pi /
                                     period * static_cast<float>(n)) +
                   0.435f * dsp::sin(5 * 2 * dsp::constants<float>::pi /
                                     period * static_cast<float>(n));
        }

        dsp::BufferContext ctxt(x, kN, buffer);
        CTXTRUN(ctxt) { nsdf.process(ctxt, dl); };

        auto peakPeriod = findPeakPeriod(nsdf);
        REQUIRE_THAT(peakPeriod, Catch::Matchers::WithinAbs(period, 5e-2));
    }
}
