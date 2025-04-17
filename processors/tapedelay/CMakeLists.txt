set(tapedelay_sources
    TapeDelay.cpp)
set(tapedelay_sc_sources
    TapeDelay_SC.cpp)

option(DSP_TAPEDELAY_SWITCH_INDICATOR
    "Compile Tapedelay with switch indicator functionality" OFF)

add_library(dsp_tapedelay_processor STATIC ${tapedelay_sources})
set_property(TARGET dsp_tapedelay_processor PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(dsp_tapedelay_processor
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(dsp_tapedelay_processor
    PUBLIC dsp_include
    PRIVATE dsp_compile_options
)
if(DSP_TAPEDELAY_SWITCH_INDICATOR)
    target_compile_definitions(dsp_tapedelay_processor
        PUBLIC TAPEDELAY_SWITCH_INDICATOR)
endif()

if(BUILD_SC)
target_sources(SCProcessors PRIVATE
    ${tapedelay_sc_sources})
target_link_libraries(SCProcessors PRIVATE dsp_tapedelay_processor)
endif()
