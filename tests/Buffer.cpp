#include "dsp/Buffer.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

TEST_CASE("Buffer Class", "[buffer]")
{
    using mf          = dsp::mfloat<2>;
    constexpr auto kN = 100;
    dsp::Buffer<mf, kN> buf;
    mf data[decltype(buf)::kSize]{};

#define loop(i) for (int i = 0; i < 2; ++i)

    SECTION("Constructor")
    {
        REQUIRE((buf.getData() == nullptr));
        buf.setData(data);
        REQUIRE((buf.getData() == data));
    }

    SECTION("Write/Read")
    {
        buf.setData(data);
        auto off = GENERATE(take(1, random(0, kN)));

        auto incr = GENERATE(1, 43, 64, 98);

        auto id = buf.getId();

        buf.nextId(off);
        id += off;

        mf x = {2, 3}, y;
        buf.write(0, x.load());

        REQUIRE((id == buf.getId()));
        buf.nextId(incr);
        y.store(buf.read(incr));

        loop(i)
        {
            REQUIRE((x[i] == y[i]));
            REQUIRE(((x[i] == data[id & decltype(buf)::kMask][i])));
        }
    }

    SECTION("Vectorized")
    {
        buf.setData(data);

        using mfv = dsp::mfloat<>;
#define loopv(i) for (size_t i = 0; i < mfv::kSize; ++i)

        mfv x, y, z;

        loopv(i) x[i] = static_cast<float>(i) + 1.f;

        auto off  = GENERATE(take(1, random(0, kN - 1)));
        auto incr = GENERATE(take(3, random(0, kN - 1)));

        auto id = buf.getId();

        buf.nextId(off);
        id += off;

        buf.writeVec(0, x.load());

        REQUIRE((id == buf.getId()));

        // delay with partial vector
        y.store(buf.readVec(1));

        // whole delay
        buf.nextId(incr);
        z.store(buf.readVec(incr));

        loopv(i)
        {
            // partial vector delay
            REQUIRE((y[i] == (i < 2 ? 0 : x[i - 2])));
            // whole delay
            REQUIRE((x[i] == z[i]));
            // check data
            REQUIRE(
                ((x[i] == data[(id & decltype(buf)::kMask) + i / 2][i % 2])));
        }
    }

    SECTION("Vectorized Limits")
    {
        buf.setData(data);

        using mfv = dsp::mfloat<>;

        mfv x, y;

        loopv(i) x[i] = static_cast<float>(i) + 1.f;

        auto id = buf.getId();
        buf.writeVec(id, x.load());

        // partial vector delay at end of buffer
        y.store(buf.readVec(id + 1));

        loopv(i)
        {
            // partial vector delay
            REQUIRE((y[i] == (i < 2 ? 0 : x[i - 2])));
        }
    }
}

TEST_CASE("Buffer Context", "[buffer][context]")
{

    using ft = float;

    constexpr auto kBufN      = 20;
    constexpr auto kN         = 50;
    constexpr auto kBlockSize = 10;

    dsp::Buffer<ft, kBufN> buffer;
    ft bufdata[decltype(buffer)::kSize]{};
    buffer.setData(bufdata);

    ft data[kN];
    for (size_t i = 0; i < kN; ++i) {
        data[i] = static_cast<float>(i);
    }

    SECTION("Construction")
    {
        {
            dsp::BufferContext ctxt(data, kBlockSize, buffer);

            REQUIRE((ctxt.getData() == data));
            REQUIRE((ctxt.getBuffer().getData() == bufdata));
            REQUIRE((ctxt.getBuffer().getId() == buffer.getId()));
            REQUIRE((ctxt.getBlockSize() == kBlockSize));
        }
        {
            dsp::Context ctxt(data, kBlockSize);
            dsp::BufferContext bctxt(ctxt, buffer);

            REQUIRE((bctxt.getData() == data));
            REQUIRE((bctxt.getBuffer().getData() == bufdata));
            REQUIRE((bctxt.getBuffer().getId() == buffer.getId()));
            REQUIRE((bctxt.getBlockSize() == kBlockSize));
        }
    }

    SECTION("Read/Write")
    {
        int i = 0;

        for (size_t j = 0; j < kN / kBlockSize; ++j) {
            dsp::BufferContext ctxt(data + kBlockSize * j, kBlockSize, buffer);

            CTXTRUN(ctxt)
            {
                ctxt.write(0, ctxt.getInput());

                auto delay = GENERATE(take(1, random(0, kBufN - 1)));

                auto v = ctxt.read(delay);

                if (i < delay) REQUIRE((v == 0.f));
                else
                    REQUIRE((v == data[i - delay]));
                ++i;
            };
            buffer.nextBlock(ctxt);
        }
    }

    SECTION("Read/Write Vectorized")
    {
        int i      = 0;
        auto delay = GENERATE(8, 13, kBufN - 1);

        for (size_t j = 0; j < kN / kBlockSize; ++j) {
            dsp::BufferContext ctxt(data + kBlockSize * j, kBlockSize, buffer);

            CTXTRUNVEC(ctxt)
            {
                auto v = ctxt.read(delay);
                ctxt.write(0, ctxt.getInput());

                if constexpr (std::is_same_v<float, decltype(v)>) {
                    if (i < delay) REQUIRE((v == 0.f));
                    else
                        REQUIRE((v == data[i - delay]));
                    ++i;
                } else {
                    for (size_t j = 0; j < decltype(v)::kSize; ++j) {

                        if (i < delay) REQUIRE((v[j] == 0.f));
                        else
                            REQUIRE((v[j] == data[i - delay]));
                        ++i;
                    }
                }
            };
            buffer.nextBlock(ctxt, true);
        }
    }
}
