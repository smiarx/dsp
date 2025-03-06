set(filters_sources
    Filters.cpp)

add_library(Filters MODULE ${filters_sources})

target_include_directories(Filters PRIVATE
    /usr/include/SuperCollider/common
    /usr/include/SuperCollider/server
    /usr/include/SuperCollider/plugin_interface
    ${CMAKE_SOURCE_DIR}/include
)

target_compile_options(Filters PRIVATE
       -Wall -Werror -fopenmp-simd
       $<$<CONFIG:RELEASE>:-Ofast -march=native -mtune=native -ffast-math>
       $<$<CONFIG:DEBUG>:-O0>
       $<$<CONFIG:DEBUG>:-g>
   )

