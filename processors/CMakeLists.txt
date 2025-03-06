
if(BUILD_SC)

add_library(SCProcessors MODULE SC.cpp)

target_include_directories(SCProcessors PRIVATE
    /usr/include/SuperCollider/common
    /usr/include/SuperCollider/server
    /usr/include/SuperCollider/plugin_interface
    ${CMAKE_SOURCE_DIR}/include
)

endif()

add_subdirectory(tapedelay/)
add_subdirectory(springs/)
add_subdirectory(filters/)
