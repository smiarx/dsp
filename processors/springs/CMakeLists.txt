set(springs_sources
    Springs.cpp)
set(springs_sc_sources
    Springs_SC.cpp)

option(DSP_SPRINGS_RMS "Compile Springs with RMS functionality" OFF)
option(DSP_SPRINGS_SHAKE "Compile Springs with Shake functionality" OFF)

# add avx2 version for dispatch
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    list(APPEND springs_sources "avx2/Springs.cpp")
    if(MSVC)
        set_source_files_properties(avx2/Springs.cpp COMPILE_FLAGS "/arch:AVX2")
    else()
        set_source_files_properties(avx2/Springs.cpp COMPILE_FLAGS "-mavx2 -mfma")
    endif()
endif()

add_library(dsp_springs_processor STATIC ${springs_sources})
set_property(TARGET dsp_springs_processor PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(dsp_springs_processor
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(dsp_springs_processor
    PUBLIC dsp_include
    PRIVATE dsp_compile_options
)
if(DSP_SPRINGS_RMS)
    target_compile_definitions(dsp_springs_processor
        PUBLIC SPRINGS_RMS)
endif()
if(DSP_SPRINGS_SHAKE)
    target_compile_definitions(dsp_springs_processor
        PUBLIC SPRINGS_SHAKE)
endif()

if(BUILD_SC)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    list(APPEND springs_sc_sources "avx2/Springs_SC.cpp")
    # add avx2 & fma to avx2 file
    if(MSVC)
        set_source_files_properties(avx2/Springs_SC.cpp DIRECTORY ../ PROPERTIES COMPILE_FLAGS "/arch:AVX2")
    else()
        set_source_files_properties(avx2/Springs_SC.cpp DIRECTORY ../ PROPERTIES COMPILE_FLAGS "-mavx2 -mfma")
    endif()
endif()

target_sources(SCProcessors PRIVATE
    ${springs_sc_sources})
target_link_libraries(SCProcessors PRIVATE dsp_springs_processor)

if (SUPERNOVA)
target_sources(SCProcessors_supernova PRIVATE
    ${springs_sc_sources})
target_link_libraries(SCProcessors_supernova PRIVATE dsp_springs_processor)
endif()

endif()
