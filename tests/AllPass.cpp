#include "dsp/AllPass.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

template <size_t N> dsp::fSample<N> rms(dsp::fSample<N> *x, size_t length)
{
    dsp::fSample<N> rms = {0.f};
    for (size_t n = 0; n < length; ++n) {
        for (size_t i = 0; i < N; ++i) {
            rms[i] += x[n][i] * x[n][i];
        }
    }
    for (size_t i = 0; i < N; ++i) {
        rms[i] = sqrtf(rms[i] / length);
    }

    return rms;
}

TEST_CASE("Allpass2", "[dsp][allpass2]")
{
    SECTION("Reference")
    {
        constexpr size_t N = 2;
        constexpr auto a0  = 0.3333333333333333;
        constexpr auto a1  = -0.9335804264972017;
        dsp::AllPass2<N> ap;
        decltype(ap)::State apstate{};
        ap.setCoeffs({a0, a0}, {a1, a1});

        const dsp::fData<1> expect[] = {
            0.3333333432674408,      -0.8298492431640625,
            -0.14408577978610992,    0.0972621738910675,
            0.16909800469875336,     0.1780681014060974,
            0.1652885526418686,      0.14639081060886383,
            0.12712723016738892,     0.10944774746894836,
            0.09386195242404938,     0.08035431802272797,
            0.06873562932014465,     0.058775559067726135,
            0.05025041103363037,     0.042958538979291916,
            0.03672352060675621,     0.03139297291636467,
            0.026835978031158447,    0.02294040471315384,
            0.019610285758972168,    0.01676357537508011,
            0.014330095611512661,    0.012249870225787163,
            0.010471615940332413,    0.008951504714787006,
            0.007652064319700003,    0.006541253998875618,
            0.00559169240295887,     0.004779973067343235,
            0.004086088854819536,    0.003492931602522731,
            0.002985881408676505,    0.002552437363192439,
            0.002181913238018751,    0.0018651759019121528,
            0.001594417612068355,    0.0013629640452563763,
            0.0011651098029688,      0.0009959766175597906,
            0.0008513956563547254,   0.0007278031553141773,
            0.0006221516523510218,   0.0005318369367159903,
            0.00045463297283276916,  0.00038863637018948793,
            0.00033221993362531066,  0.00028399325674399734,
            0.0002427674480713904,   0.00020752614364027977,
            0.00017740059411153197,  0.00015164828801061958,
            0.00012963431072421372,  0.00011081594857387245,
            0.00009472934470977634,  0.00008097795216599479,
            0.00006922279862919822,  0.00005917408998357132,
            0.00005058410170022398,  0.000043241059756837785,
            0.00003696398562169634,  0.00003159810512443073,
            0.000027011161364498548,
        };

        dsp::fSample<N> x{};

        for (size_t n = 0; n < sizeof(expect) / sizeof(expect[0]); ++n) {
            x[0] = n == 0;
            x[1] = n == 0;
            ap.process(dsp::Context(&x, 1), apstate);
            REQUIRE_THAT(x[0], Catch::Matchers::WithinAbs(expect[n][0], 1e-6));
            REQUIRE_THAT(x[1], Catch::Matchers::WithinAbs(expect[n][0], 1e-6));
        }
    }

    SECTION("Phasing")
    {
        constexpr size_t N           = 2;
        constexpr auto sR            = 48000.f;
        constexpr size_t PrepareSize = 0.5 * sR;
        constexpr size_t TestSize    = 64;
        float freq = GENERATE(take(2, random(200.f, sR * 0.9f)));
        float R    = GENERATE(take(2, random(0.1f, 10.f)));
        float f    = 2.f * freq / sR;

        dsp::AllPass2<N> ap;
        decltype(ap)::State apstate{};
        ap.setFreq({f, f}, {R, R});

        dsp::fSample<N> x;
        for (size_t n = 0; n < PrepareSize; ++n) {
            for (size_t i = 0; i < N; ++i) {
                x[i] = sin(M_PIf * f * n);
            }
            ap.process(dsp::Context(&x, 1), apstate);
        }

        for (size_t n = PrepareSize; n < PrepareSize + TestSize; ++n) {
            for (size_t i = 0; i < N; ++i) {
                x[i] = sin(M_PIf * f * n);
            }
            auto x0 = x;
            ap.process(dsp::Context(&x, 1), apstate);

            arrayFor(x, i)
            {
                // output is phased at 180°
                x0[i] += x[i];
                REQUIRE_THAT(x0[i], Catch::Matchers::WithinAbs(0.f, 4e-3));
            }
        }
    }

    SECTION("Unchanged power")
    {
        constexpr size_t N    = 2;
        constexpr auto sR     = 48000.f;
        constexpr size_t Size = 2048;
        float sinefreq        = GENERATE(take(2, random(200.f, sR * 0.5f)));
        float freq            = GENERATE(take(2, random(200.f, sR * 0.5f)));
        float R               = GENERATE(take(2, random(0.1f, 10.f)));
        float f               = 2.f * freq / sR;
        float sinef           = 2.f * sinefreq / sR;

        dsp::AllPass2<N> ap;
        decltype(ap)::State apstate{};
        ap.setFreq({f, f}, {R, R});

        dsp::fSample<N> x[Size];
        for (size_t n = 0; n < Size; ++n) {
            for (size_t i = 0; i < N; ++i) {
                x[n][i] = sin(M_PIf * f * n);
            }
        }

        auto sinRMS = rms(&x[Size / 2], Size / 2);
        ap.processBlock(dsp::Context(x, Size), apstate);
        auto apRMS = rms(&x[Size / 2], Size / 2);

        for (size_t i = 0; i < N; ++i) {
            REQUIRE_THAT(apRMS[i], Catch::Matchers::WithinAbs(sinRMS[i], 1e-3));
        }
    }

    //    SECTION("Energy Preserving")
    //    {
    //        constexpr size_t N    = 2;
    //        constexpr auto sR     = 48000.f;
    //        constexpr size_t Size = 1<<15;
    //        float sinefreq            = GENERATE(take(2, random(200.f, sR *
    //        0.5f))); float freq            = GENERATE(take(2, random(200.f, sR
    //        * 0.5f))); float R               = GENERATE(take(2,
    //        random(0.1f, 10.f))); float f               = 2.f * freq / sR;
    //        float of               = 1.f - f;
    //        float oR               = R+0.8f;
    //
    //        dsp::AllPass2<N,true> ap;
    //        decltype(ap)::State apstate{};
    //        ap.setFreq({f, f}, {R, R});
    //
    //        dsp::fSample<N> x[Size];
    //        for (size_t n = 0; n < Size; ++n) {
    //            for (size_t i = 0; i < N; ++i) {
    //                x[n][i] = sin(M_PIf * sinefreq * n);
    //            }
    //        }
    //
    //        auto sinRMS = rms(&x[Size / 2], Size / 2);
    //
    //        contextFor(dsp::Context(x,3*Size/4))
    //        {
    //            //ap.pregain(c);
    //            ap.process(c, apstate);
    //            //ap.postgain(c);
    //        }
    //
    //        ap.setFreq({of, of}, {oR, oR});
    //        contextFor(dsp::Context(x+3*Size/4,Size/4))
    //        {
    //            //ap.pregain(c);
    //            ap.process(c, apstate);
    //            //ap.postgain(c);
    //        }
    //
    //        auto apRMS = rms(&x[Size / 2], Size / 2);
    //
    //        for (size_t i = 0; i < N; ++i) {
    //            REQUIRE_THAT(apRMS[i], Catch::Matchers::WithinAbs(sinRMS[i],
    //            2e-2));
    //        }
    //    }
}
