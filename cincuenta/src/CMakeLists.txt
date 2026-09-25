# This is repeated here but should possibly just use the custom command and info
# from dmrgpp/src

add_executable(cincuenta cincuenta.cpp)
target_link_libraries(cincuenta PUBLIC dmrgpp_utils lanczos_globals dmrg_runner kronutil)

add_executable(freeGimp freeGimp.cpp)
target_link_libraries(freeGimp PUBLIC dmrgpp_utils kronutil)

install(TARGETS cincuenta freeGimp RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Build the unit tests Missing for old work still minimal for new dev
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

# Executable Level Tests
set(PATH_TO_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/inputs)

# exact impurity solver
add_test(NAME cincuenta0 COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputCincuenta0.ain) # energy 0.456977
# dmrg impurity solver, 3-bath fast CI check — energy -4.86046
# Keep the run isolated from all MPI variants because cincuenta uses fixed names
# for correction-vector and restart artifacts.
set(CINCUENTA_DMRG_1_DIR ${CMAKE_CURRENT_BINARY_DIR}/cincuentaDmrg3BathRank1)
file(MAKE_DIRECTORY ${CINCUENTA_DMRG_1_DIR})

if(DMRGPP_USE_MPI)
  # Launch even a one-rank MPI build through mpiexec.  Besides exercising the
  # same initialization path as multi-rank jobs, this avoids Open MPI singleton
  # teardown problems in some Kokkos-enabled environments.
  set(_cincuenta_rank1_cmd "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1")
  list(APPEND _cincuenta_rank1_cmd ${MPIEXEC_PREFLAGS})
  list(APPEND _cincuenta_rank1_cmd $<TARGET_FILE:cincuenta>)
  list(APPEND _cincuenta_rank1_cmd ${MPIEXEC_POSTFLAGS})
  list(APPEND _cincuenta_rank1_cmd -I ${PATH_TO_INPUTS} -f inputCincuentaDmrg3Bath.ain)
  add_test(NAME cincuentaDmrg3Bath COMMAND ${_cincuenta_rank1_cmd})
  set_tests_properties(cincuentaDmrg3Bath PROPERTIES LABELS "MPI" PROCESSORS 1 WORKING_DIRECTORY
                                                     ${CINCUENTA_DMRG_1_DIR})
else()
  add_test(NAME cincuentaDmrg3Bath COMMAND $<TARGET_FILE:cincuenta> -I ${PATH_TO_INPUTS} -f inputCincuentaDmrg3Bath.ain)
  set_tests_properties(cincuentaDmrg3Bath PROPERTIES WORKING_DIRECTORY ${CINCUENTA_DMRG_1_DIR})
endif()

if(DMRGPP_USE_MPI)
  set(CINCUENTA_DMRG_MPI_DIR ${CMAKE_CURRENT_BINARY_DIR}/cincuentaDmrg3BathMpi)
  file(MAKE_DIRECTORY ${CINCUENTA_DMRG_MPI_DIR})

  set(_cincuenta_mpi_cmd "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2")
  list(APPEND _cincuenta_mpi_cmd ${MPIEXEC_PREFLAGS})
  list(APPEND _cincuenta_mpi_cmd $<TARGET_FILE:cincuenta>)
  list(APPEND _cincuenta_mpi_cmd ${MPIEXEC_POSTFLAGS})
  list(APPEND _cincuenta_mpi_cmd -I ${PATH_TO_INPUTS} -f inputCincuentaDmrg3Bath.ain)
  add_test(NAME cincuentaDmrg3Bath_mpi COMMAND ${_cincuenta_mpi_cmd})
  set_tests_properties(cincuentaDmrg3Bath_mpi PROPERTIES LABELS "MPI" PROCESSORS 2 WORKING_DIRECTORY
                                                         ${CINCUENTA_DMRG_MPI_DIR})

  # Check a Matsubara value assembled on rank zero after both distributed
  # particle/hole runs.  A successful check also proves the shared artifact was written.
  add_check(cincuentaDmrg3Bath_mpi_gimp "0\\.0628319 [-0-9.eE+]+ ([-0-9.eE+]+)" "-0.58646" "1e-4"
            ${CINCUENTA_DMRG_MPI_DIR}/gimp_dmrg.txt cincuentaDmrg3Bath_mpi)
  set_tests_properties(cincuentaDmrg3Bath_mpi_gimp PROPERTIES LABELS "MPI")

  # Interacting equilibrium DMFT preflight: two complete self-consistency
  # iterations, followed by the final real-frequency solve.  Rank counts use
  # separate directories, and their assembled numerical outputs are compared.
  if(Python3_FOUND)
    set(CINCUENTA_SCF2_RANK1_DIR ${CMAKE_CURRENT_BINARY_DIR}/cincuentaDmrg3BathScf2Rank1)
    set(CINCUENTA_SCF2_RANK2_DIR ${CMAKE_CURRENT_BINARY_DIR}/cincuentaDmrg3BathScf2Rank2)
    file(MAKE_DIRECTORY ${CINCUENTA_SCF2_RANK1_DIR} ${CINCUENTA_SCF2_RANK2_DIR})

    set(_cincuenta_scf2_rank1_cmd "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1")
    list(APPEND _cincuenta_scf2_rank1_cmd ${MPIEXEC_PREFLAGS})
    list(APPEND _cincuenta_scf2_rank1_cmd $<TARGET_FILE:cincuenta>)
    list(APPEND _cincuenta_scf2_rank1_cmd ${MPIEXEC_POSTFLAGS})
    list(APPEND _cincuenta_scf2_rank1_cmd -I ${PATH_TO_INPUTS} -f inputCincuentaDmrg3BathScf2.ain)
    add_test(NAME cincuentaDmrg3BathScf2_rank1 COMMAND ${_cincuenta_scf2_rank1_cmd})
    set_tests_properties(
      cincuentaDmrg3BathScf2_rank1
      PROPERTIES FIXTURES_SETUP
                 CincuentaScf2Rank1
                 LABELS
                 "MPI;DMFTPreflight"
                 PROCESSORS
                 1
                 TIMEOUT
                 300
                 WORKING_DIRECTORY
                 ${CINCUENTA_SCF2_RANK1_DIR})

    set(_cincuenta_scf2_rank2_cmd "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2")
    list(APPEND _cincuenta_scf2_rank2_cmd ${MPIEXEC_PREFLAGS})
    list(APPEND _cincuenta_scf2_rank2_cmd $<TARGET_FILE:cincuenta>)
    list(APPEND _cincuenta_scf2_rank2_cmd ${MPIEXEC_POSTFLAGS})
    list(APPEND _cincuenta_scf2_rank2_cmd -I ${PATH_TO_INPUTS} -f inputCincuentaDmrg3BathScf2.ain)
    add_test(NAME cincuentaDmrg3BathScf2_rank2 COMMAND ${_cincuenta_scf2_rank2_cmd})
    set_tests_properties(
      cincuentaDmrg3BathScf2_rank2
      PROPERTIES FIXTURES_SETUP
                 CincuentaScf2Rank2
                 LABELS
                 "MPI;DMFTPreflight"
                 PROCESSORS
                 2
                 TIMEOUT
                 300
                 WORKING_DIRECTORY
                 ${CINCUENTA_SCF2_RANK2_DIR})

    set(_scf2_log runForinputCincuentaDmrg3BathScf2.cout)
    add_check(cincuentaDmrg3BathScf2_rank1_iterations "SelfConsistLoop iter= ([0-9]+)" "1" "0"
              ${CINCUENTA_SCF2_RANK1_DIR}/${_scf2_log} cincuentaDmrg3BathScf2_rank1)
    add_check(cincuentaDmrg3BathScf2_rank2_iterations "SelfConsistLoop iter= ([0-9]+)" "1" "0"
              ${CINCUENTA_SCF2_RANK2_DIR}/${_scf2_log} cincuentaDmrg3BathScf2_rank2)
    set_tests_properties(cincuentaDmrg3BathScf2_rank1_iterations cincuentaDmrg3BathScf2_rank2_iterations
                         PROPERTIES LABELS "MPI;DMFTPreflight")

    add_test(
      NAME cincuentaDmrg3BathScf2_compare_rank2
      COMMAND
        ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/compare_dmft_mpi_outputs.py
        ${CINCUENTA_SCF2_RANK1_DIR} ${CINCUENTA_SCF2_RANK2_DIR} --tolerance 1e-5 --log ${_scf2_log} --iterations 2)
    set_tests_properties(
      cincuentaDmrg3BathScf2_compare_rank2 PROPERTIES FIXTURES_REQUIRED "CincuentaScf2Rank1;CincuentaScf2Rank2" LABELS
                                                      "MPI;DMFTPreflight")

    option(CINCUENTA_MPI_PREFLIGHT_4 "Enable the four-rank cincuenta DMFT preflight" OFF)
    if(CINCUENTA_MPI_PREFLIGHT_4)
      set(CINCUENTA_SCF2_RANK4_DIR ${CMAKE_CURRENT_BINARY_DIR}/cincuentaDmrg3BathScf2Rank4)
      file(MAKE_DIRECTORY ${CINCUENTA_SCF2_RANK4_DIR})
      set(_cincuenta_scf2_rank4_cmd "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4")
      list(APPEND _cincuenta_scf2_rank4_cmd ${MPIEXEC_PREFLAGS})
      list(APPEND _cincuenta_scf2_rank4_cmd $<TARGET_FILE:cincuenta>)
      list(APPEND _cincuenta_scf2_rank4_cmd ${MPIEXEC_POSTFLAGS})
      list(APPEND _cincuenta_scf2_rank4_cmd -I ${PATH_TO_INPUTS} -f inputCincuentaDmrg3BathScf2.ain)
      add_test(NAME cincuentaDmrg3BathScf2_rank4 COMMAND ${_cincuenta_scf2_rank4_cmd})
      set_tests_properties(
        cincuentaDmrg3BathScf2_rank4
        PROPERTIES FIXTURES_SETUP
                   CincuentaScf2Rank4
                   LABELS
                   "MPI;DMFTPreflight"
                   PROCESSORS
                   4
                   TIMEOUT
                   300
                   WORKING_DIRECTORY
                   ${CINCUENTA_SCF2_RANK4_DIR})
      add_check(cincuentaDmrg3BathScf2_rank4_iterations "SelfConsistLoop iter= ([0-9]+)" "1" "0"
                ${CINCUENTA_SCF2_RANK4_DIR}/${_scf2_log} cincuentaDmrg3BathScf2_rank4)
      set_tests_properties(cincuentaDmrg3BathScf2_rank4_iterations PROPERTIES LABELS "MPI;DMFTPreflight")
      add_test(
        NAME cincuentaDmrg3BathScf2_compare_rank4
        COMMAND
          ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/compare_dmft_mpi_outputs.py
          ${CINCUENTA_SCF2_RANK1_DIR} ${CINCUENTA_SCF2_RANK4_DIR} --tolerance 1e-5 --log ${_scf2_log} --iterations 2)
      set_tests_properties(
        cincuentaDmrg3BathScf2_compare_rank4 PROPERTIES FIXTURES_REQUIRED "CincuentaScf2Rank1;CincuentaScf2Rank4"
                                                        LABELS "MPI;DMFTPreflight")
    endif()
  endif()
endif()

# dmrg impurity solver — slow; excluded from fast CI, run with: ctest -L Nightly
add_test(NAME cincuenta1 COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputCincuenta1.ain) # energy 0.456977
set_tests_properties(cincuenta1 PROPERTIES LABELS "Nightly")

# Large fitted-first-bath/rank-3 stress case. GBEK paper-figure and report
# generation starts from the true atomic limit and must not depend on this run.
set(NEQ_GBEK_NEAR_ATOMIC_DIR ${CMAKE_CURRENT_BINARY_DIR}/neqGbekFig3L3NearAtomic)
file(MAKE_DIRECTORY ${NEQ_GBEK_NEAR_ATOMIC_DIR})
add_test(
  NAME neqGBEKFig3L3NearAtomic
  COMMAND $<TARGET_FILE:cincuenta> -I ${PATH_TO_INPUTS} -f inputNeqGBEKFig3L3.ain
  WORKING_DIRECTORY ${NEQ_GBEK_NEAR_ATOMIC_DIR})
set_tests_properties(neqGBEKFig3L3NearAtomic PROPERTIES LABELS "Nightly" RUN_SERIAL TRUE)

# ---- U=0 self-consistency tests: G_imp must equal G_cluster_0 ----------------
# Physical basis: at U=0 the self-energy Sigma=0 exactly, so the DMFT condition
# G_imp = local lattice GF must hold for any correct impurity solver. cincuenta
# writes gimp_<solver>.txt and latticeG_<solver>.txt after the Matsubara loop;
# compare_gimp_latticeG.py checks they agree to within bath-fitting error.
if(Python3_FOUND)
  if(NOT PATH_TO_TESTSUITE)
    set(PATH_TO_TESTSUITE ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite)
  endif()

  # -- ExactDiag ---------------------------------------------------------------
  set(U0_ED_DIR ${CMAKE_CURRENT_BINARY_DIR}/u0ExactDiag)
  file(MAKE_DIRECTORY ${U0_ED_DIR})

  add_test(
    NAME u0_gImpEqGcluster_exactdiag
    COMMAND $<TARGET_FILE:cincuenta> -I ${PATH_TO_INPUTS} -f inputU0ExactDiag.ain
    WORKING_DIRECTORY ${U0_ED_DIR})

  add_test(NAME u0_gImpEqGcluster_exactdiag_check
           COMMAND ${Python3_EXECUTABLE} ${PATH_TO_TESTSUITE}/compare_gimp_latticeG.py ${U0_ED_DIR}/gimp_exactdiag.txt
                   ${U0_ED_DIR}/latticeG_exactdiag.txt --tolerance 0.015)
  set_tests_properties(u0_gImpEqGcluster_exactdiag_check PROPERTIES DEPENDS u0_gImpEqGcluster_exactdiag)

  # -- DMRG --------------------------------------------------------------------
  set(U0_DMRG_DIR ${CMAKE_CURRENT_BINARY_DIR}/u0Dmrg)
  file(MAKE_DIRECTORY ${U0_DMRG_DIR})

  add_test(
    NAME u0_gImpEqGcluster_dmrg
    COMMAND $<TARGET_FILE:cincuenta> -I ${PATH_TO_INPUTS} -f inputU0Dmrg.ain
    WORKING_DIRECTORY ${U0_DMRG_DIR})

  add_test(NAME u0_gImpEqGcluster_dmrg_check
           COMMAND ${Python3_EXECUTABLE} ${PATH_TO_TESTSUITE}/compare_gimp_latticeG.py ${U0_DMRG_DIR}/gimp_dmrg.txt
                   ${U0_DMRG_DIR}/latticeG_dmrg.txt --tolerance 0.05)
  set_tests_properties(u0_gImpEqGcluster_dmrg_check PROPERTIES DEPENDS u0_gImpEqGcluster_dmrg LABELS "Nightly")
  set_tests_properties(u0_gImpEqGcluster_dmrg PROPERTIES LABELS "Nightly")

endif()

# Fast end-to-end GBEK regression: a fitted non-atomic first bath plus a
# positive-rank (L=2) Cholesky second bath.  Keep the output isolated so this
# test cannot overwrite other non-equilibrium runs.
if(Python3_FOUND)
  set(NEQ_GBEK_FAST_DIR ${CMAKE_CURRENT_BINARY_DIR}/neqGbekFast)
  file(MAKE_DIRECTORY ${NEQ_GBEK_FAST_DIR})

  add_test(
    NAME neqGBEKFast
    COMMAND $<TARGET_FILE:cincuenta> -I ${PATH_TO_INPUTS} -f inputNeqGBEKFast.ain
    WORKING_DIRECTORY ${NEQ_GBEK_FAST_DIR})
  set_tests_properties(neqGBEKFast PROPERTIES FIXTURES_SETUP neq_gbek_fast)

  add_test(NAME neqGBEKFast_check
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/check_neq_gbek_fast.py
                   ${NEQ_GBEK_FAST_DIR} --matsubaras 16)
  set_tests_properties(neqGBEKFast_check PROPERTIES FIXTURES_REQUIRED neq_gbek_fast)
endif()
