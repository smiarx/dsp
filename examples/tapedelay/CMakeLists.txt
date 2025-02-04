cmake_minimum_required(VERSION 3.15)

project(Test VERSION 0.0.1)

set(sources
    TapeDelay.cpp)

add_library(SCTapeDelay MODULE ${sources})

target_include_directories(SCTapeDelay PRIVATE
    /usr/include/SuperCollider/common
    /usr/include/SuperCollider/server
    /usr/include/SuperCollider/plugin_interface
)

target_compile_options(SCTapeDelay PRIVATE
       -Wall -Werror -fopenmp-simd
       $<$<CONFIG:RELEASE>:-Ofast -g -march=native -mtune=native -ffast-math>
       $<$<CONFIG:DEBUG>:-Og>
       $<$<CONFIG:DEBUG>:-g>
   )

