
if(BUILD_SC)

set(SC_TARGET SCProcessors)

add_library(${SC_TARGET} MODULE SC.cpp)

target_include_directories(SCProcessors PRIVATE
    /usr/include/SuperCollider/common
    /usr/include/SuperCollider/server
    /usr/include/SuperCollider/plugin_interface
    ${CMAKE_SOURCE_DIR}/include
)

target_compile_options(${SC_TARGET} PRIVATE -Wall -Wextra -Wpedantic -Werror -fopenmp-simd)

endif()

add_subdirectory(tapedelay/)
add_subdirectory(springs/)
add_subdirectory(filters/)
