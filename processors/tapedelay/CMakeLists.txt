option(DSP_TAPEDELAY_SWITCH_INDICATOR
    "Compile Tapedelay with switch indicator functionality" OFF)

set(tapedelay_sources
    TapeDelay.cpp)

#add sc file to source if compile for SC
if(BUILD_SC)
    list(APPEND tapedelay_sources "TapeDelay_SC.cpp")
endif()

# main library
add_library(dsp_tapedelay_processor INTERFACE)

add_library(dsp_tapedelay_base INTERFACE)
target_include_directories(dsp_tapedelay_base
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(dsp_tapedelay_base
    INTERFACE
    dsp_include
)
if(BUILD_SC)
    target_link_libraries(dsp_tapedelay_base
        INTERFACE SuperCollider)
endif()
if(DSP_TAPEDELAY_SWITCH_INDICATOR)
    target_compile_definitions(dsp_tapedelay_base
        INTERFACE TAPEDELAY_SWITCH_INDICATOR)
endif()

add_library(dsp_tapedelay_generic STATIC ${tapedelay_sources})
set_property(TARGET dsp_tapedelay_generic PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET dsp_tapedelay_generic PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
target_link_libraries(dsp_tapedelay_generic
    PUBLIC  dsp_tapedelay_base
    PRIVATE dsp_compile_options
)
target_link_libraries(dsp_tapedelay_processor
    INTERFACE dsp_tapedelay_generic)

if(DSP_X86_DISPATCH)
### AVX2 library
add_library(dsp_tapedelay_avx2 STATIC ${tapedelay_sources})
set_property(TARGET dsp_tapedelay_avx2 PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET dsp_tapedelay_avx2 PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
target_link_libraries(dsp_tapedelay_avx2
    PUBLIC  dsp_tapedelay_base
    PRIVATE dsp_compile_options dsp_avx2
)
target_link_libraries(dsp_tapedelay_processor
    INTERFACE dsp_tapedelay_avx2)
endif()

if(BUILD_SC)
    target_link_libraries(${SC_TARGET}
        PUBLIC dsp_tapedelay_processor)

    if (SUPERNOVA)
        target_link_libraries(${SC_TARGET}_supernova
            PUBLIC dsp_tapedelay_processor)
    endif()
endif()
