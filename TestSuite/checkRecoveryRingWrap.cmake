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
file(READ "${INPUT_FILE}" input_text)
string(REPLACE "OutputFile=\"data210\";" "OutputFile=\"data210ring\";" ring_input "${input_text}")
string(REPLACE "@M=200" "@M=2" ring_input "${ring_input}")
set(ring_input_file "${WORK_ROOT}/input210-ring.ain")
file(WRITE "${ring_input_file}" "${ring_input}")

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${ring_input_file}"
  WORKING_DIRECTORY "${WORK_ROOT}"
  RESULT_VARIABLE initial_result
  OUTPUT_VARIABLE initial_stdout
  ERROR_VARIABLE initial_stderr)
if(NOT initial_result EQUAL 0)
  message(FATAL_ERROR "Initial ring-wrap RecoverySave run failed:\n${initial_stdout}${initial_stderr}")
endif()

set(latest_file "${WORK_ROOT}/Recovery0data210ring.hd5")
set(previous_file "${WORK_ROOT}/Recovery1data210ring.hd5")
foreach(recovery_file IN ITEMS "${latest_file}" "${previous_file}")
  if(NOT EXISTS "${recovery_file}")
    message(FATAL_ERROR "Expected ring-wrap recovery file is missing: ${recovery_file}")
  endif()
endforeach()

file(SHA256 "${WORK_ROOT}/data210ring.hd5" output_before)
file(SHA256 "${latest_file}" latest_before)
file(SHA256 "${previous_file}" previous_before)

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${ring_input_file}"
  WORKING_DIRECTORY "${WORK_ROOT}"
  RESULT_VARIABLE rerun_result
  OUTPUT_VARIABLE rerun_stdout
  ERROR_VARIABLE rerun_stderr)
set(rerun_output "${rerun_stdout}${rerun_stderr}")
if(NOT rerun_result EQUAL 0)
  message(FATAL_ERROR "Ring-wrap RecoverySave rerun failed:\n${rerun_output}")
endif()
if(NOT rerun_output MATCHES "RecoverySave: checkpoint already completed; exiting\\.")
  message(FATAL_ERROR "Ring-wrap RecoverySave rerun did not select the terminal generation:\n${rerun_output}")
endif()

file(SHA256 "${WORK_ROOT}/data210ring.hd5" output_after)
file(SHA256 "${latest_file}" latest_after)
file(SHA256 "${previous_file}" previous_after)
if(NOT output_before STREQUAL output_after
   OR NOT latest_before STREQUAL latest_after
   OR NOT previous_before STREQUAL previous_after)
  message(FATAL_ERROR "Ring-wrap RecoverySave rerun modified a completed artifact")
endif()

message(STATUS "Ring-wrap RecoverySave selects the most recent generation")
