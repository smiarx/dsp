cmake_minimum_required(VERSION 3.15)

project(Test VERSION 0.0.1)

set(sources
    Reverb.cpp)

add_library(SCReverb MODULE ${sources})

target_include_directories(SCReverb PRIVATE
    /usr/include/SuperCollider/common
    /usr/include/SuperCollider/server
    /usr/include/SuperCollider/plugin_interface
)
