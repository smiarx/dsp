cmake_minimum_required(VERSION 3.15)

project(dsp VERSION 0.0.1)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    add_subdirectory(tests)
    add_subdirectory(bench)
endif()

OPTION(BUILD_SC  "Build SuperCollider processors"  ON)
add_subdirectory(processors/)
