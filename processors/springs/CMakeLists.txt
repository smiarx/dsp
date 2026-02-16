option(DSP_SPRINGS_RMS "Compile Springs with RMS functionality" OFF)
option(DSP_SPRINGS_SHAKE "Compile Springs with Shake functionality" OFF)

set(springs_sources
    Springs.cpp)

#add sc file to source if compile for SC
if(BUILD_SC)
    list(APPEND springs_sources "Springs_SC.cpp")
endif()

# main library
add_library(dsp_springs_processor INTERFACE)

add_library(dsp_springs_base INTERFACE)
target_include_directories(dsp_springs_base
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(dsp_springs_base
    INTERFACE
    dsp_include
)
if(BUILD_SC)
    target_link_libraries(dsp_springs_base
        INTERFACE SuperCollider)
endif()
if(DSP_SPRINGS_RMS)
    target_compile_definitions(dsp_springs_base
        INTERFACE SPRINGS_RMS)
endif()
if(DSP_SPRINGS_SHAKE)
    target_compile_definitions(dsp_springs_base
        INTERFACE SPRINGS_SHAKE)
endif()

add_library(dsp_springs_generic STATIC ${springs_sources})
set_property(TARGET dsp_springs_generic PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET dsp_springs_generic PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
target_link_libraries(dsp_springs_generic
    PUBLIC  dsp_springs_base
    PRIVATE dsp_compile_options
)
target_link_libraries(dsp_springs_processor
    INTERFACE dsp_springs_generic)

if(DSP_X86_DISPATCH)
### AVX2 library
add_library(dsp_springs_avx2 STATIC ${springs_sources})
set_property(TARGET dsp_springs_avx2 PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET dsp_springs_avx2 PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
target_link_libraries(dsp_springs_avx2
    PUBLIC  dsp_springs_base
    PRIVATE dsp_compile_options dsp_avx2
)
target_link_libraries(dsp_springs_processor
    INTERFACE dsp_springs_avx2)
endif()

if(BUILD_SC)
    target_link_libraries(${SC_TARGET}
        PUBLIC dsp_springs_processor)

    if (SUPERNOVA)
        target_link_libraries(${SC_TARGET}_supernova
            PUBLIC dsp_springs_processor)
    endif()
endif()
