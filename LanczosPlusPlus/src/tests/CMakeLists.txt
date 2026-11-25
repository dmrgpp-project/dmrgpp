add_executable(test_Combinatorial test_Combinatorial.cpp)
target_include_directories(test_Combinatorial PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(test_Combinatorial PRIVATE Catch2::Catch2WithMain lanczos_globals)
catch_discover_tests(test_Combinatorial)

add_executable(test_BasisOneSpin test_BasisOneSpin.cpp)
target_include_directories(test_BasisOneSpin PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(test_BasisOneSpin PRIVATE Catch2::Catch2WithMain lanczos_globals)
catch_discover_tests(test_BasisOneSpin)

add_executable(test_EngineSpectral test_EngineSpectral.cpp)
target_include_directories(test_EngineSpectral PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(test_EngineSpectral PRIVATE Catch2::Catch2 lanczos_globals)
# Engine construction uses MPI in MPI-enabled builds. Initialize it in the
# test main and launch the suite through mpiexec even for one rank.
if(DMRGPP_USE_MPI)
  add_test(
    NAME "Engine spectral regression tests"
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:test_EngineSpectral> ${MPIEXEC_POSTFLAGS})
  set_tests_properties("Engine spectral regression tests" PROPERTIES LABELS "MPI;EngineSpectral" PROCESSORS 1)
else()
  catch_discover_tests(test_EngineSpectral)
endif()
