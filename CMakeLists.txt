cmake_minimum_required(VERSION 3.25)

project(dsp VERSION 0.0.1)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    OPTION(BUILD_SC  "Build SuperCollider processors"  ON)
    option(BUILD_TESTING "Build tests" ON)
else()
    set(BUILD_SC OFF)
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

add_library(dsp_compile_options INTERFACE)
if(MSVC)
    target_compile_options(dsp_compile_options INTERFACE
    /openmp
    $<$<CONFIG:RELEASE>: /fp:fast>
)
else()
    target_compile_options(dsp_compile_options INTERFACE
    -Wall -Wextra -Wpedantic -Werror -fopenmp-simd
    $<$<CONFIG:RELEASE>:-O3 -ffast-math -funsafe-math-optimizations -ffp-contract=fast>
)
endif()
target_compile_features(dsp_compile_options INTERFACE cxx_std_17)

add_subdirectory(processors/)
