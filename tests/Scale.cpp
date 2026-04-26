#include <catch2/catch_all.hpp>
#include <dsp/music/Scale.hpp>

using namespace dsp; // Assuming Scale is under 'dsp' namespace

TEST_CASE("Scale", "[music][scale]")
{

    using enum Note;
    static constexpr auto minor = A | B | C | D | E | F | G;
    static constexpr auto major = A | B | Cs | D | E | Fs | Gs;
    SECTION("Operations")
    {
        static_assert(major.contains(A));
        static_assert(major.contains(B));
        static_assert(major.contains(Cs));
        static_assert(major.contains(D));
        static_assert(major.contains(E));
        static_assert(major.contains(Fs));
        static_assert(major.contains(Gs));
        static_assert(!major.contains(As));
        static_assert(!major.contains(C));
        static_assert(!major.contains(Ds));
        static_assert(!major.contains(F));
        static_assert(!major.contains(G));

        static_assert(minor != major);
        constexpr auto intersect = minor & major;
        static_assert(intersect.contains(A));
        static_assert(intersect.contains(B));
        static_assert(intersect.contains(D));
        static_assert(intersect.contains(E));
        static_assert(!intersect.contains(As));
        static_assert(!intersect.contains(C));
        static_assert(!intersect.contains(Cs));
        static_assert(!intersect.contains(Ds));
        static_assert(!intersect.contains(F));
        static_assert(!intersect.contains(Fs));
        static_assert(!intersect.contains(Gs));
        static_assert(!intersect.contains(G));

        Note note   = GENERATE(A, B, C, D, E, F, G);
        auto scale2 = major ^ note;
        REQUIRE(major.contains(note) != scale2.contains(note));

        static_assert((~minor) == (As | Cs | Ds | Fs | Gs));
    }

    SECTION("Transpose")
    {
        constexpr auto transposed = minor >> A << E;
        constexpr auto minorE     = E | Fs | G | A | B | C | D;
        static_assert(transposed == minorE);
    }

    SECTION("forEachNote")
    {
        minor.forEachNote([](Note note) {
            REQUIRE(minor.contains(note));
            REQUIRE(note != As);
            REQUIRE(note != Cs);
            REQUIRE(note != Ds);
            REQUIRE(note != Fs);
            REQUIRE(note != Gs);
        });
    }
}
