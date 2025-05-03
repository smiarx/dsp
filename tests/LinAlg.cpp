#include "dsp/LinAlg.h"
#include "dsp/Signal.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>

TEST_CASE("Linear Algebra", "[dsp][linalg]")
{
    SECTION("Matrix mult vector")
    {
        constexpr size_t kN              = 4;
        dsp::linalg::Matrix<float, kN> m = {{{
            {0.9606811539349068, 0.28769739052079435, 0.43830855665372603,
             0.07081992794185832},

            {0.3671050259076819, 0.03221069684247624, 0.48051192620455274,
             0.877288993134758},

            {0.29311101955536323, 0.4163668941981783, 0.5987133508019171,
             0.6853979608797864},

            {0.9613021771529163, 0.9649997209773631, 0.12888976383945794,
             0.5058441074502963},
        }}};

        dsp::Sample<float, kN> x      = {0.20873822, 0.08766247, 0.47080132,
                                         0.14881192};
        dsp::Sample<float, kN> expect = {0.51376248, 0.40250666, 0.43466998,
                                         0.48965005};

        auto result = m.mult(x);

        REQUIRE_THAT(result[0], Catch::Matchers::WithinAbs(expect[0], 1e-6));
        REQUIRE_THAT(result[1], Catch::Matchers::WithinAbs(expect[1], 1e-6));
        REQUIRE_THAT(result[2], Catch::Matchers::WithinAbs(expect[2], 1e-6));
        REQUIRE_THAT(result[3], Catch::Matchers::WithinAbs(expect[3], 1e-6));
    }

    SECTION("Matrix mult vector")
    {
        constexpr size_t kN               = 4;
        dsp::linalg::Matrix<float, kN> m1 = {{{
            {0.0910283384099071, 0.48987322566042635, 0.3477358717973138,
             0.14085546647215907},

            {0.05965172967910248, 0.3662754999849357, 0.2258265068658456,
             0.39373537797612357},

            {0.53364789051963, 0.05103040062499764, 0.008064225605900477,
             0.3766549704122426},

            {0.6937804176459286, 0.4811429022446442, 0.5321093552441208,
             0.8997876496657035},
        }}};

        dsp::linalg::Matrix<float, kN> m2 = {{{
            {0.41708431543756763, 0.09330762039197038, 0.7590547759006182,
             0.8952000403914135},

            {0.72722868657717, 0.5181388010248434, 0.18017720769473555,
             0.872264239292175},

            {0.5062309229672677, 0.6859722360078181, 0.4087377229526331,
             0.31669567803917065},

            {0.8546789935652717, 0.6052328259522658, 0.8221567169377872,
             0.9331338128296616},
        }}};

        dsp::linalg::Matrix<float, kN> expect = {{{
            {1.06967269100736, 0.7079487491212193, 0.6485720172860362,
             1.1868788114748163},

            {0.7984173297178734, 0.9749096736027327, 0.8354859284541799,
             1.159160943138745},

            {0.5248400735652778, 0.6724777262955022, 0.5027582514549355,
             0.7803088851700942},

            {1.2000353576396172, 1.1312920089443301, 0.9370394486096871,
             1.507979477974491},
        }}};

        auto result = m1.mult(m2);

        REQUIRE_THAT(result[0][0],
                     Catch::Matchers::WithinAbs(expect[0][0], 1e-6));
        REQUIRE_THAT(result[0][1],
                     Catch::Matchers::WithinAbs(expect[0][1], 1e-6));
        REQUIRE_THAT(result[0][2],
                     Catch::Matchers::WithinAbs(expect[0][2], 1e-6));
        REQUIRE_THAT(result[0][3],
                     Catch::Matchers::WithinAbs(expect[0][3], 1e-6));

        REQUIRE_THAT(result[1][0],
                     Catch::Matchers::WithinAbs(expect[1][0], 1e-6));
        REQUIRE_THAT(result[1][1],
                     Catch::Matchers::WithinAbs(expect[1][1], 1e-6));
        REQUIRE_THAT(result[1][2],
                     Catch::Matchers::WithinAbs(expect[1][2], 1e-6));
        REQUIRE_THAT(result[1][3],
                     Catch::Matchers::WithinAbs(expect[1][3], 1e-6));

        REQUIRE_THAT(result[2][0],
                     Catch::Matchers::WithinAbs(expect[2][0], 1e-6));
        REQUIRE_THAT(result[2][1],
                     Catch::Matchers::WithinAbs(expect[2][1], 1e-6));
        REQUIRE_THAT(result[2][2],
                     Catch::Matchers::WithinAbs(expect[2][2], 1e-6));
        REQUIRE_THAT(result[2][3],
                     Catch::Matchers::WithinAbs(expect[2][3], 1e-6));

        REQUIRE_THAT(result[3][0],
                     Catch::Matchers::WithinAbs(expect[3][0], 1e-6));
        REQUIRE_THAT(result[3][1],
                     Catch::Matchers::WithinAbs(expect[3][1], 1e-6));
        REQUIRE_THAT(result[3][2],
                     Catch::Matchers::WithinAbs(expect[3][2], 1e-6));
        REQUIRE_THAT(result[3][3],
                     Catch::Matchers::WithinAbs(expect[3][3], 1e-6));
    }

    SECTION("identity")
    {
        constexpr auto kN = 4;
        auto id           = dsp::linalg::identity<float, kN>();
        auto x            = dsp::fSample<kN>{};
        x[0]              = GENERATE(take(2, random(0, 1)));
        x[1]              = GENERATE(take(2, random(0, 1)));
        x[2]              = GENERATE(take(2, random(0, 1)));
        x[3]              = GENERATE(take(2, random(0, 1)));

        auto r = id.mult(x);

        REQUIRE_THAT(r[0], Catch::Matchers::WithinAbs(x[0], 1e-6));
        REQUIRE_THAT(r[1], Catch::Matchers::WithinAbs(x[1], 1e-6));
        REQUIRE_THAT(r[2], Catch::Matchers::WithinAbs(x[2], 1e-6));
        REQUIRE_THAT(r[3], Catch::Matchers::WithinAbs(x[3], 1e-6));
    }
}
