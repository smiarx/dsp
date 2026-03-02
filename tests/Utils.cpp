#include <catch2/catch_test_macros.hpp>
#include <dsp/utils/Utils.hpp>

TEST_CASE("Modular Arithmetic")
{
    constexpr auto kMod = 256;
    REQUIRE(dsp::modCmpGT<kMod>(3, 0));
    REQUIRE(dsp::modCmpGT<kMod>(3, 255));
    REQUIRE(dsp::modCmpGT<kMod>(100, 250));
    REQUIRE(dsp::modCmpGT<kMod>(20, 263));
    REQUIRE(!dsp::modCmpGT<kMod>(210, 20));
    REQUIRE(!dsp::modCmpGT<kMod>(255, 0));
    REQUIRE(dsp::modCmpGT<kMod>(10, 9));
    REQUIRE(!dsp::modCmpGT<kMod>(10, 10));
    REQUIRE(!dsp::modCmpGT<kMod>(10, 11));
    REQUIRE(dsp::modCmpGE<kMod>(256, 0));
    REQUIRE(dsp::modCmpGE<kMod>(356, 90));
    REQUIRE(dsp::modCmpGE<kMod>(30, 260));
    REQUIRE(dsp::modCmpGE<kMod>(20, 19));
    REQUIRE(dsp::modCmpGE<kMod>(20, 20));
    REQUIRE(!dsp::modCmpGE<kMod>(20, 21));
    REQUIRE(dsp::modCmpEq<kMod>(0, 256));
    REQUIRE(dsp::modCmpEq<kMod>(60, 316));
    REQUIRE((dsp::modDiff<kMod>(3, 250) == 9));
    REQUIRE((dsp::modDiff<kMod>(240, 8) == -24));
    REQUIRE((dsp::modAdd<kMod>(240, 7) == 247));
    REQUIRE((dsp::modAdd<kMod>(4, -17) == 243));
}
