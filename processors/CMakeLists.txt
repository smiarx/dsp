if(BUILD_SC)

set(SC_TARGET SCProcessors)

add_library(${SC_TARGET} MODULE SC.cpp)

add_library(SuperCollider INTERFACE)
target_include_directories(SuperCollider INTERFACE
    /usr/include/SuperCollider/common
    /usr/include/SuperCollider/server
    /usr/include/SuperCollider/plugin_interface
)
set_property(TARGET SuperCollider PROPERTY SYSTEM TRUE)

target_link_libraries(${SC_TARGET}
    PUBLIC dsp_include
    PRIVATE dsp_compile_options SuperCollider
)

endif()

add_subdirectory(tapedelay/)
#add_subdirectory(springs/)
add_subdirectory(filters/)
