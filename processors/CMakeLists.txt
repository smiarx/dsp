if(BUILD_SC)

option(SUPERNOVA "Build plugins for supernova" OFF)

set(SC_TARGET SCProcessors)

add_library(${SC_TARGET} MODULE SC.cpp)
target_link_libraries(${SC_TARGET}
    PUBLIC dsp_include
    PRIVATE dsp_compile_options
)

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

if (SUPERNOVA)
add_library(${SC_TARGET}_supernova MODULE SC.cpp)
target_link_libraries(${SC_TARGET}_supernova
    PUBLIC dsp_include
    PRIVATE dsp_compile_options
)

target_compile_definitions(${SC_TARGET}_supernova
    PUBLIC SUPERNOVA)

target_include_directories(SuperCollider INTERFACE /usr/include/SuperCollider/external_libraries/nova-tt)
# actually just boost.atomic
target_include_directories(SuperCollider INTERFACE /usr/include/SuperCollider/external_libraries/boost)
target_include_directories(SuperCollider INTERFACE /usr/include/SuperCollider/external_libraries/boost_lockfree)
target_include_directories(SuperCollider INTERFACE /usr/include/SuperCollider/external_libraries/boost-lockfree)

target_link_libraries(${SC_TARGET}_supernova
    PUBLIC dsp_include
    PRIVATE dsp_compile_options SuperCollider
)
endif()

endif()

add_subdirectory(tapedelay/)
add_subdirectory(springs/)
add_subdirectory(voicebox/)
add_subdirectory(filters/)
# need modifications
#add_subdirectory(adaptive/)
