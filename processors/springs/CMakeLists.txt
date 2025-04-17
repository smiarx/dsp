set(springs_sources
    Springs.cpp)
set(springs_sc_sources
    Springs_SC.cpp)

option(DSP_SPRINGS_RMS "Compile Springs with RMS functionality" OFF)

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

if(BUILD_SC)
target_sources(SCProcessors PRIVATE
    ${springs_sc_sources})
target_link_libraries(SCProcessors PRIVATE dsp_springs_processor)
endif()
