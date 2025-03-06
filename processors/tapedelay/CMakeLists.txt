set(tapedelay_sources
    TapeDelay.cpp)
set(tapedelay_sc_sources
    TapeDelay_SC.cpp)

if(BUILD_SC)
target_sources(SCProcessors PRIVATE
    ${tapedelay_sources}
    ${tapedelay_sc_sources})
endif()
