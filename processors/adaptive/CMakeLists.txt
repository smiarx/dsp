set(adaptive_sc_sources
    Adaptive_SC.cpp)

if(BUILD_SC)
target_sources(SCProcessors PRIVATE
    ${adaptive_sc_sources})

if (SUPERNOVA)
target_sources(SCProcessors_supernova PRIVATE
    ${adaptive_sc_sources})
endif()

endif()
