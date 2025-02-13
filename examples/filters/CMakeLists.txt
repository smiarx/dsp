cmake_minimum_required(VERSION 3.15)

project(Filters VERSION 0.0.1)

set(sources
    Filters.cpp)

add_library(Filters MODULE ${sources})

target_include_directories(Filters PRIVATE
    /usr/include/SuperCollider/common
    /usr/include/SuperCollider/server
    /usr/include/SuperCollider/plugin_interface
)

target_compile_options(Filters PRIVATE
       -Wall -Werror -fopenmp-simd
       $<$<CONFIG:RELEASE>:-Ofast -march=native -mtune=native -ffast-math>
       $<$<CONFIG:DEBUG>:-O0>
       $<$<CONFIG:DEBUG>:-g>
   )

