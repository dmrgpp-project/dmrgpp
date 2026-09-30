add_executable(test_Fit test_Fit.cpp)
target_link_libraries(test_Fit PRIVATE Catch2::Catch2WithMain dmrgpp_utils psimaglite)
target_include_directories(test_Fit PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(test_Fit)

add_executable(test_NeqBathDecomposition test_NeqBathDecomposition.cpp)
target_link_libraries(test_NeqBathDecomposition PRIVATE Catch2::Catch2WithMain psimaglite)
target_include_directories(test_NeqBathDecomposition PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src)
catch_discover_tests(test_NeqBathDecomposition)

add_executable(test_NeqRealTimeGf test_NeqRealTimeGf.cpp)
target_link_libraries(test_NeqRealTimeGf PRIVATE Catch2::Catch2WithMain psimaglite)
target_include_directories(test_NeqRealTimeGf PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src)
catch_discover_tests(test_NeqRealTimeGf)

add_executable(test_spectralWeight test_spectralWeight.cpp)
target_link_libraries(test_spectralWeight PRIVATE Catch2::Catch2WithMain psimaglite dmrgpp_utils)
target_include_directories(test_spectralWeight PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src
                                                       ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(test_spectralWeight)

add_executable(test_ImpuritySolverNeqGBEK test_ImpuritySolverNeqGBEK.cpp)
target_link_libraries(test_ImpuritySolverNeqGBEK PRIVATE Catch2::Catch2WithMain dmrgpp_utils lanczos_globals
                                                         dmrg_runner kronutil)
target_include_directories(
  test_ImpuritySolverNeqGBEK PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src ${CMAKE_SOURCE_DIR}/dmrg/Engine
                                     ${CMAKE_SOURCE_DIR}/LanczosPlusPlus/src/Engine ${CMAKE_SOURCE_DIR})
catch_discover_tests(test_ImpuritySolverNeqGBEK)

add_executable(test_NeqGbekInput test_NeqGbekInput.cpp)
target_link_libraries(test_NeqGbekInput PRIVATE Catch2::Catch2WithMain dmrgpp_utils lanczos_globals dmrg_runner
                                                kronutil)
target_include_directories(test_NeqGbekInput PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src ${CMAKE_SOURCE_DIR}/dmrg/Engine
                                                     ${CMAKE_SOURCE_DIR}/LanczosPlusPlus/src/Engine ${CMAKE_SOURCE_DIR})
catch_discover_tests(test_NeqGbekInput)

add_executable(test_DmftSolverEquilibriumInitialData test_DmftSolverEquilibriumInitialData.cpp)
target_link_libraries(test_DmftSolverEquilibriumInitialData PRIVATE Catch2::Catch2WithMain dmrgpp_utils lanczos_globals
                                                                    dmrg_runner kronutil)
target_include_directories(
  test_DmftSolverEquilibriumInitialData PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src ${CMAKE_SOURCE_DIR}/dmrg/Engine
                                                ${CMAKE_SOURCE_DIR}/LanczosPlusPlus/src/Engine ${CMAKE_SOURCE_DIR})
set(DMFT_EQUILIBRIUM_INITIAL_DATA_DIR ${CMAKE_CURRENT_BINARY_DIR}/equilibriumInitialData)
file(MAKE_DIRECTORY ${DMFT_EQUILIBRIUM_INITIAL_DATA_DIR})
# Launch this MPI-linked test through mpiexec even at one rank.  Running it
# as an Open MPI singleton can abort during teardown in Kokkos-enabled builds.
if(DMRGPP_USE_MPI)
  set(_equilibrium_initial_data_cmd "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1")
  list(APPEND _equilibrium_initial_data_cmd ${MPIEXEC_PREFLAGS})
  list(APPEND _equilibrium_initial_data_cmd $<TARGET_FILE:test_DmftSolverEquilibriumInitialData>)
  list(APPEND _equilibrium_initial_data_cmd ${MPIEXEC_POSTFLAGS})
  add_test(NAME test_DmftSolverEquilibriumInitialData COMMAND ${_equilibrium_initial_data_cmd})
  set_tests_properties(test_DmftSolverEquilibriumInitialData PROPERTIES PROCESSORS 1)
else()
  add_test(NAME test_DmftSolverEquilibriumInitialData COMMAND test_DmftSolverEquilibriumInitialData)
endif()
set_tests_properties(test_DmftSolverEquilibriumInitialData PROPERTIES WORKING_DIRECTORY ${DMFT_EQUILIBRIUM_INITIAL_DATA_DIR})
