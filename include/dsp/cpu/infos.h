#pragma once

#include "defines.h"
#include <cstdint>

#if defined(_MSC_VER)
// Contains the definition of __cpuidex
#include <intrin.h>
#endif

namespace dsp::cpu
{

#if DSP_X86_64

namespace internal
{
static auto xcr0() noexcept
{
    uint32_t xcr0;
#if defined(_MSC_VER) && _MSC_VER >= 1400
    xcr0 = (uint32_t)_xgetbv(0);
#elif defined(__GNUC__)
    __asm__("xorl %%ecx, %%ecx\n"
            "xgetbv\n"
            : "=a"(xcr0)
            :
#if defined(__i386__)
            : "ecx", "edx"
#else
            : "rcx", "rdx"
#endif
    );
#else /* _MSC_VER < 1400 */
#error "_MSC_VER < 1400 is not supported"
#endif /* _MSC_VER && _MSC_VER >= 1400 */
    return xcr0;
};

static auto cpuid(int reg[4], unsigned int level,
                  unsigned int count = 0) noexcept
{
#if defined(_MSC_VER)
    __cpuidex(reg, level, count);
#elif defined(__INTEL_COMPILER)
    __cpuid(reg, level);
#elif defined(__GNUC__) || defined(__clang__)
#if defined(__i386__) && defined(__PIC__)
    // %ebx may be the PIC register
    __asm__("xchg{l}\t{%%}ebx, %1\n\t"
            "cpuid\n\t"
            "xchg{l}\t{%%}ebx, %1\n\t"
            : "=a"(reg[0]), "=r"(reg[1]), "=c"(reg[2]), "=d"(reg[3])
            : "0"(level), "2"(count));

#else
    __asm__("cpuid\n\t"
            : "=a"(reg[0]), "=b"(reg[1]), "=c"(reg[2]), "=d"(reg[3])
            : "0"(level), "2"(count));
#endif
#else
#error "Unsupported configuration"
#endif
};

} // namespace internal

struct Infos {
    uint32_t sse2 : 1;
    uint32_t sse3 : 1;
    uint32_t ssse3 : 1;
    uint32_t sse4_1 : 1;
    uint32_t sse4_2 : 1;
    uint32_t fma3_sse42 : 1;
    uint32_t avx : 1;
    uint32_t fma4 : 1;
    uint32_t avx2 : 1;
    uint32_t avxvnni : 1;
};

static inline auto getInfos() noexcept
{
    int regs1[4];

    internal::cpuid(regs1, 0x1);

    // OS can explicitly disable the usage of SSE/AVX extensions
    // by setting an appropriate flag in CR0 register
    //
    // https://docs.kernel.org/admin-guide/hw-vuln/gather_data_sampling.html

    unsigned sseStateOsEnabled = 1;
    unsigned avxStateOsEnabled = 1;

    // OSXSAVE: A value of 1 indicates that the OS has set CR4.OSXSAVE[bit
    // 18] to enable XSETBV/XGETBV instructions to access XCR0 and
    // to support processor extended state management using
    // XSAVE/XRSTOR.
    bool osxsave = regs1[2] >> 27 & 1;
    if (osxsave) {
        uint32_t xcr0 = internal::xcr0();

        sseStateOsEnabled = xcr0 >> 1 & 1;
        avxStateOsEnabled = xcr0 >> 2 & sseStateOsEnabled;
    }

    Infos infos;

    infos.sse2       = regs1[3] >> 26 & sseStateOsEnabled;
    infos.sse3       = regs1[2] >> 0 & sseStateOsEnabled;
    infos.ssse3      = regs1[2] >> 9 & sseStateOsEnabled;
    infos.sse4_1     = regs1[2] >> 19 & sseStateOsEnabled;
    infos.sse4_2     = regs1[2] >> 20 & sseStateOsEnabled;
    infos.fma3_sse42 = regs1[2] >> 12 & sseStateOsEnabled;

    infos.avx = regs1[2] >> 28 & avxStateOsEnabled;

    int regs8[4];
    internal::cpuid(regs8, 0x80000001);
    infos.fma4 = regs8[2] >> 16 & avxStateOsEnabled;

    int regs7[4];
    internal::cpuid(regs7, 0x7);
    infos.avx2 = regs7[1] >> 5 & avxStateOsEnabled;

    int regs7a[4];
    internal::cpuid(regs7a, 0x7, 0x1);
    infos.avxvnni = regs7a[0] >> 4 & avxStateOsEnabled;

    return infos;
}

#endif
} // namespace dsp::cpu
