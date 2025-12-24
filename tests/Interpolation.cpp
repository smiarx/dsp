#include "dsp/Interpolation.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("Interpolation")
{
    using namespace dsp::interp;
    using namespace Catch::Matchers;

    SECTION("Polynomial")
    {
        double y, e;
        SECTION("2nd degree")
        {
            SECTION("")
            {
                y = Polynomial{0., 1., 1.}(1.);
                e = 2.;
            }
            SECTION("")
            {
                y = Polynomial{0.473, 1.2389, -2.382}(0.398);
                e = 0.5887638719999999;
            }
        }
        SECTION("3rd degree")
        {
            y = Polynomial{1.382, -0.38623, 2.38372, -0.28372}(1.2932);
            e = 4.255377832613527;
        }
        REQUIRE_THAT(y, WithinAbs(e, 1e-10));
    }
    SECTION("Parabolic")
    {
        auto poly = LagrangeParabolic{1., 2., 3.};
        REQUIRE_THAT(poly(-1.), WithinAbs(1., 1e-10));
        REQUIRE_THAT(poly(0.), WithinAbs(2., 1e-10));
        REQUIRE_THAT(poly(1.), WithinAbs(3., 1e-10));
    }
    SECTION("Cubic")
    {
        auto poly = LagrangeCubic{1., 2., 3., 4.};
        REQUIRE_THAT(poly(-1.), WithinAbs(1., 1e-10));
        REQUIRE_THAT(poly(0.), WithinAbs(2., 1e-10));
        REQUIRE_THAT(poly(1.), WithinAbs(3., 1e-10));
        REQUIRE_THAT(poly(2.), WithinAbs(4., 1e-10));
    }
    SECTION("Hermite")
    {
        auto poly = Hermite{1., 2., 3., 4.};
        REQUIRE_THAT(poly(0.), WithinAbs(2., 1e-10));
        REQUIRE_THAT(poly(1.), WithinAbs(3., 1e-10));

        // check tangents
        auto dx = 0.0001;

        auto dp0  = (poly(dx) - poly(-dx)) / (2 * dx);
        auto edp0 = (3. - 1.) / 2.;
        REQUIRE_THAT(dp0, WithinAbs(edp0, 1e-10));

        auto dp1  = (poly(1. + dx) - poly(1. - dx)) / (2 * dx);
        auto edp1 = (4. - 2.) / 2.;
        REQUIRE_THAT(dp1, WithinAbs(edp1, 1e-10));
    }

    SECTION("Parabolic extremum")
    {
        auto ym1 = 0.288706, y0 = 0.925619, y1 = 0.91472;
        auto ex = 0.48317567442406123, ey = 1.0012376841622261;

        auto [maxx, maxy] = extrema(LagrangeParabolic(ym1, y0, y1));
        REQUIRE_THAT(maxx, WithinAbs(ex, 1e-9));
        REQUIRE_THAT(maxy, WithinAbs(ey, 1e-9));
    }

    SECTION("Cubic extrema")
    {
        double ym1, y0, y1, y2, expectx, expecty;

        ym1 = 0.455381, y0 = 0.931866, y1 = 0.933195, y2 = 0.473251;
        expectx = 0.49913210865052404, expecty = 0.9910574888469426;

        auto [p1, p2]       = extrema(Hermite{ym1, y0, y1, y2});
        auto [extrx, extry] = p1;
        REQUIRE_THAT(extrx, WithinAbs(expectx, 1e-9));
        REQUIRE_THAT(extry, WithinAbs(expecty, 1e-9));
    }
}
