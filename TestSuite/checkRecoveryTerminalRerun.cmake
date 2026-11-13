if(NOT DEFINED DMRG_EXECUTABLE)
  message(FATAL_ERROR "DMRG_EXECUTABLE is required")
endif()

if(NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE is required")
endif()

if(NOT DEFINED WORK_ROOT)
  message(FATAL_ERROR "WORK_ROOT is required")
endif()

file(REMOVE_RECURSE "${WORK_ROOT}")
file(MAKE_DIRECTORY "${WORK_ROOT}")

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${INPUT_FILE}"
  WORKING_DIRECTORY "${WORK_ROOT}"
  RESULT_VARIABLE initial_result
  OUTPUT_VARIABLE initial_stdout
  ERROR_VARIABLE initial_stderr)
if(NOT initial_result EQUAL 0)
  message(FATAL_ERROR "Initial retained RecoverySave run failed:\n${initial_stdout}${initial_stderr}")
endif()

set(output_file "${WORK_ROOT}/data210.hd5")
if(NOT EXISTS "${output_file}")
  message(FATAL_ERROR "Completed output is missing: ${output_file}")
endif()
file(SHA256 "${output_file}" output_before)
file(GLOB recovery_files "${WORK_ROOT}/Recovery*data210.hd5")
if(NOT recovery_files)
  message(FATAL_ERROR "No retained RecoverySave files found")
endif()
set(recovery_before)
foreach(recovery_file IN LISTS recovery_files)
  file(SHA256 "${recovery_file}" recovery_hash)
  list(APPEND recovery_before "${recovery_file}:${recovery_hash}")
endforeach()

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${INPUT_FILE}"
  WORKING_DIRECTORY "${WORK_ROOT}"
  RESULT_VARIABLE rerun_result
  OUTPUT_VARIABLE rerun_stdout
  ERROR_VARIABLE rerun_stderr)
set(rerun_output "${rerun_stdout}${rerun_stderr}")
if(NOT rerun_result EQUAL 0)
  message(FATAL_ERROR "Terminal RecoverySave rerun failed:\n${rerun_output}")
endif()
if(NOT rerun_output MATCHES "RecoverySave: checkpoint already completed; exiting\\.")
  message(FATAL_ERROR "Terminal RecoverySave rerun did not report completion:\n${rerun_output}")
endif()

file(SHA256 "${output_file}" output_after)
if(NOT output_before STREQUAL output_after)
  message(FATAL_ERROR "Terminal RecoverySave rerun modified ${output_file}")
endif()
foreach(recovery_entry IN LISTS recovery_before)
  string(REPLACE ":" ";" recovery_parts "${recovery_entry}")
  list(GET recovery_parts 0 recovery_file)
  list(GET recovery_parts 1 expected_hash)
  if(NOT EXISTS "${recovery_file}")
    message(FATAL_ERROR "Terminal RecoverySave rerun removed ${recovery_file}")
  endif()
  file(SHA256 "${recovery_file}" actual_hash)
  if(NOT expected_hash STREQUAL actual_hash)
    message(FATAL_ERROR "Terminal RecoverySave rerun modified ${recovery_file}")
  endif()
endforeach()

message(STATUS "Terminal RecoverySave rerun is an idempotent no-op")
