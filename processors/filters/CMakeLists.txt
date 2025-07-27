set(filters_sc_sources
    Filters_SC.cpp)

if(BUILD_SC)
target_sources(SCProcessors PRIVATE
    ${filters_sc_sources})

if (SUPERNOVA)
target_sources(SCProcessors_supernova PRIVATE
    ${filters_sc_sources})
endif()

endif()
