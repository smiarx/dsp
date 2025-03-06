set(springs_sources
    Springs.cpp)
set(springs_sc_sources
    Springs_SC.cpp)

if(BUILD_SC)
target_sources(SCProcessors PRIVATE
    ${springs_sources}
    ${springs_sc_sources})
endif()
