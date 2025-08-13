set(voicebox_sources
    VoiceBox.cpp)
set(voicebox_sc_sources
    VoiceBox_SC.cpp)

add_library(dsp_voicebox_processor STATIC ${voicebox_sources})
set_property(TARGET dsp_voicebox_processor PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(dsp_voicebox_processor
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(dsp_voicebox_processor
    PUBLIC dsp_include
    PRIVATE dsp_compile_options
)

if(BUILD_SC)
target_sources(SCProcessors PRIVATE
    ${voicebox_sc_sources})
target_link_libraries(SCProcessors PRIVATE dsp_voicebox_processor)

if (SUPERNOVA)
target_sources(SCProcessors_supernova PRIVATE
    ${voicebox_sc_sources})
target_link_libraries(SCProcessors_supernova PRIVATE dsp_voicebox_processor)
endif()

endif()
