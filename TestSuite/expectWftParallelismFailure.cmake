if(NOT DEFINED DMRG_EXECUTABLE)
  message(FATAL_ERROR "DMRG_EXECUTABLE is required")
endif()

if(NOT DEFINED INPUT_GS)
  message(FATAL_ERROR "INPUT_GS is required")
endif()

if(NOT DEFINED INPUT_WFT)
  message(FATAL_ERROR "INPUT_WFT is required")
endif()

if(NOT DEFINED MODE)
  message(FATAL_ERROR "MODE is required")
endif()

if(NOT DEFINED WORK_DIR)
  message(FATAL_ERROR "WORK_DIR is required")
endif()

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
file(COPY "${INPUT_GS}" DESTINATION "${WORK_DIR}")
file(COPY "${INPUT_WFT}" DESTINATION "${WORK_DIR}")

get_filename_component(input_gs_name "${INPUT_GS}" NAME)
get_filename_component(input_wft_name "${INPUT_WFT}" NAME)
set(input_gs_path "${WORK_DIR}/${input_gs_name}")
set(input_wft_path "${WORK_DIR}/${input_wft_name}")

file(READ "${input_wft_path}" input_wft)
if(MODE STREQUAL "blocks")
  string(APPEND input_wft "\nThreads=2;\n")
  set(expected_diagnostic "ParallelWftMany=2 conflicts with accelerated WFT blocks section WftAccelBlocks=2")
elseif(MODE STREQUAL "noaccel")
  string(REPLACE "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart;"
                 "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart,wftNoAccel;" input_wft
                 "${input_wft}")
  string(APPEND input_wft "\nThreads=1;\nThreads:ParallelWftMany=2;\nThreads:WaveFunctionTransfLocal=2;\n")
elseif(MODE STREQUAL "patches-ctor")
  string(REPLACE "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart;"
                 "SolverOptions=geometryallinsystem,TargetingExpression,restart,wftAccelPatches;" input_wft
                 "${input_wft}")
  string(APPEND input_wft "\nThreads=1;\nThreads:ParallelWftMany=2;\nThreads:BlockDiagWfCtor=2;\n")
  set(expected_diagnostic "ParallelWftMany=2 conflicts with accelerated WFT patches section BlockDiagWfCtor=2")
elseif(MODE STREQUAL "patches-transform")
  string(REPLACE "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart;"
                 "SolverOptions=geometryallinsystem,TargetingExpression,restart,wftAccelPatches;" input_wft
                 "${input_wft}")
  string(APPEND input_wft "\nThreads=1;\nThreads:ParallelWftMany=2;\nThreads:BlockDiagWfTransform=2;\n")
  set(expected_diagnostic "ParallelWftMany=2 conflicts with accelerated WFT patches section BlockDiagWfTransform=2")
else()
  message(FATAL_ERROR "Unknown MODE=${MODE}")
endif()
file(WRITE "${input_wft_path}" "${input_wft}")

if(MODE MATCHES "^patches-")
  file(READ "${input_gs_path}" input_gs)
  string(REPLACE "SolverOptions=twositedmrg,geometryallinsystem;" "SolverOptions=geometryallinsystem;" input_gs
                 "${input_gs}")
  file(WRITE "${input_gs_path}" "${input_gs}")
endif()

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${input_gs_path}"
  WORKING_DIRECTORY "${WORK_DIR}"
  RESULT_VARIABLE gs_result
  OUTPUT_VARIABLE gs_stdout
  ERROR_VARIABLE gs_stderr)
if(NOT gs_result EQUAL 0)
  message(FATAL_ERROR "Ground-state producer failed:\n${gs_stdout}${gs_stderr}")
endif()

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${input_wft_path}"
  WORKING_DIRECTORY "${WORK_DIR}"
  RESULT_VARIABLE wft_result
  OUTPUT_VARIABLE wft_stdout
  ERROR_VARIABLE wft_stderr)
set(output "${wft_stdout}${wft_stderr}")
if(MODE STREQUAL "noaccel")
  if(NOT wft_result EQUAL 0)
    message(FATAL_ERROR "No-acceleration WFT run failed:\n${output}")
  endif()
  message(STATUS "No-acceleration WFT completed with outer and local threads")
  return()
endif()

if(wft_result EQUAL 0)
  message(FATAL_ERROR "Conflicting WFT run unexpectedly succeeded:\n${output}")
endif()

string(FIND "${output}" "${expected_diagnostic}" diagnostic_position)
if(diagnostic_position EQUAL -1)
  message(FATAL_ERROR "Expected diagnostic not found: ${expected_diagnostic}\n${output}")
endif()

message(STATUS "Observed expected failure: ${expected_diagnostic}")
