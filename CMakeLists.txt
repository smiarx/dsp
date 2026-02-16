cmake_minimum_required(VERSION 3.25)

project(dsp VERSION 0.0.1)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    OPTION(BUILD_SC  "Build SuperCollider processors"  ON)
    option(BUILD_TESTING "Build tests" ON)
else()
    set(BUILD_SC OFF)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)"
   OR APPLE)
   option(DSP_X86_DISPATCH "Build for different x86 instructions set" OFF)

   if(DSP_X86_DISPATCH)
       add_library(dsp_avx2 INTERFACE)
       if(MSVC)
           target_compile_options(dsp_avx2
               INTERFACE /arch:AVX2)
       else()
           target_compile_options(dsp_avx2
               INTERFACE -mavx2 -mfma)
       endif()
       target_compile_definitions(dsp_avx2
           INTERFACE DSP_COMPILE_AVX2=1)
   endif()
else()
    set(DSP_X86_DISPATCH OFF)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    include(CTest)
    add_subdirectory(bench)
endif()
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

add_library(dsp_include INTERFACE)
target_include_directories(dsp_include INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
#add gcem library (directly because CMAKE version isn't up to date)
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/submodules/gcem/include/>
)
if(DSP_X86_DISPATCH)
    target_compile_definitions(dsp_include
        INTERFACE
        DSP_X86_DISPATCH=1
    )
endif()

add_library(dsp_compile_options INTERFACE)
if(MSVC)
    target_compile_options(dsp_compile_options INTERFACE
    /openmp
    $<$<CONFIG:RELEASE>: /fp:fast>
)
else()
    target_compile_options(dsp_compile_options INTERFACE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wno-shadow
    $<$<CONFIG:RELEASE>: -ffast-math -funsafe-math-optimizations -ffp-contract=fast>
)
endif()
target_compile_features(dsp_compile_options INTERFACE cxx_std_17)

add_subdirectory(processors/)
