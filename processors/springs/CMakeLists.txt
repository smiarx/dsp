set(springs_sources
    Springs.cpp)

add_library(SCSprings MODULE ${springs_sources})

target_include_directories(SCSprings PRIVATE
    /usr/include/SuperCollider/common
    /usr/include/SuperCollider/server
    /usr/include/SuperCollider/plugin_interface
    ${CMAKE_SOURCE_DIR}/include
)

target_compile_options(SCSprings PRIVATE
       -Wall -Werror -fopenmp-simd
       $<$<CONFIG:RELEASE>:-Ofast -march=native -mtune=native -ffast-math>
       $<$<CONFIG:DEBUG>:-O0>
       $<$<CONFIG:DEBUG>:-g>
   )

