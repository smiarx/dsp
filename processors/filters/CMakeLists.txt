set(filters_sc_sources
    Filters_SC.cpp)

if(BUILD_SC)
target_sources(SCProcessors PRIVATE
    ${filters_sc_sources})
endif()
