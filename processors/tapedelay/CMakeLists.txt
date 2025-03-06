set(tapedelay_sources
    TapeDelay.cpp)

add_library(SCTapeDelay MODULE ${tapedelay_sources})

target_include_directories(SCTapeDelay PRIVATE
    /usr/include/SuperCollider/common
    /usr/include/SuperCollider/server
    /usr/include/SuperCollider/plugin_interface
    ${CMAKE_SOURCE_DIR}/include
)

target_compile_options(SCTapeDelay PRIVATE
       -Wall -Werror -fopenmp-simd
       $<$<CONFIG:RELEASE>:-Ofast -march=native -mtune=native -ffast-math -funroll-loops>
       $<$<CONFIG:DEBUG>:-O0>         
       $<$<CONFIG:DEBUG>:-g>
   )

