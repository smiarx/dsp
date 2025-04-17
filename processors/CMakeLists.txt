if(BUILD_SC)

set(SC_TARGET SCProcessors)

add_library(${SC_TARGET} MODULE SC.cpp)

target_include_directories(${SC_TARGET} PRIVATE
    /usr/include/SuperCollider/common
    /usr/include/SuperCollider/server
    /usr/include/SuperCollider/plugin_interface
)
target_link_libraries(${SC_TARGET}
    PUBLIC dsp_include
    PRIVATE dsp_compile_options
)

endif()

add_subdirectory(tapedelay/)
add_subdirectory(springs/)
add_subdirectory(filters/)
