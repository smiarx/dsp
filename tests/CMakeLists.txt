find_package(Catch2 3 REQUIRED)

set(test_sources
    SIMD.cpp
    Signal.cpp
    Context.cpp
    FIRFilter.cpp
)

add_executable(tests ${test_sources})
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)

