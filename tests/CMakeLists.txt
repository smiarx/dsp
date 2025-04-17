find_package(Catch2 3 REQUIRED)

set(test_sources
    SIMD.cpp
    Signal.cpp
    Context.cpp
    FIRFilter.cpp
    VAFilters.cpp
    RMS.cpp
    LinAlg.cpp
    Hadamard.cpp
    AllPass.cpp
)

add_executable(tests ${test_sources})
target_link_libraries(tests
    PUBLIC dsp_include
    PRIVATE
        dsp_compile_options
        Catch2::Catch2WithMain
)

enable_testing()
add_test(NAME tests COMMAND $<TARGET_FILE:tests>)

