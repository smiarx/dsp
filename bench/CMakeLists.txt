find_package(Catch2 3 REQUIRED)

set(bench_sources
    AllPass.cpp
    #FIRFilter.cpp
    #VAFilters.cpp
    #TapeDelay.cpp
    #../processors/tapedelay/TapeDelay.cpp
    #Springs.cpp
    #../processors/springs/Springs.cpp
)

add_executable(bench ${bench_sources})
target_include_directories(bench PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(bench PRIVATE Catch2::Catch2WithMain)
