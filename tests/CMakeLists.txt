find_package(Catch2 3 QUIET)

if(NOT ${Catch2_FOUND})
    Include(FetchContent)

    FetchContent_Declare(
      Catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
      GIT_TAG        v3.4.0 # or a later release
    )
    FetchContent_MakeAvailable(Catch2)
endif()

set(test_sources
    simd.cpp
    #Signal.cpp
    #Context.cpp
    #FIRFilter.cpp
    #VAFilters.cpp
    #RMS.cpp
    #LinAlg.cpp
    #Hadamard.cpp
    #AllPass.cpp
    #Enveloppe.cpp
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

