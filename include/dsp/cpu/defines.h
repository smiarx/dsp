#pragma once

#if defined(__x86_64__)
#define DSP_X86_64 1
#endif

#if defined(__i386__)
#define DSP_I386 1
#endif

#if defined(__aarch64__)
#define DSP_AARCH64 1
#endif

#if defined(__arm__)
#define DSP_ARM32 1
#endif

#if DSP_X86_64 || DSP_I386

#define DSP_X86 1

#if defined(__AVX2__)
#define DSP_AVX          1
#define DSP_SSE2         1
#define DSP_MAX_VEC_SIZE 32
#elif defined(__SSE2__)
#define DSP_SSE2         1
#define DSP_MAX_VEC_SIZE 16

#ifndef DSP_NO_DISPATCH
#define DSP_X86_DISPATCH 1
#endif

#endif
#endif

#if DSP_AARCH64 || DSP_ARM32
#define DSP_ARM 1

#if defined(__ARM_NEON)
#define DSP_NEON 1
#endif
#endif
