if(NOT DEFINED DMRG_EXECUTABLE)
  message(FATAL_ERROR "DMRG_EXECUTABLE is required")
endif()

if(NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE is required")
endif()

if(NOT DEFINED EXPECTED_DIAGNOSTIC)
  message(FATAL_ERROR "EXPECTED_DIAGNOSTIC is required")
endif()

set(command "${DMRG_EXECUTABLE}" -f "${INPUT_FILE}")
if(DEFINED EXTRA_OPTION)
  list(APPEND command -o "${EXTRA_OPTION}")
endif()

execute_process(
  COMMAND ${command}
  RESULT_VARIABLE result
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE stderr)

set(output "${stdout}${stderr}")
if(result EQUAL 0)
  message(FATAL_ERROR "Command unexpectedly succeeded:\n${output}")
endif()

string(FIND "${output}" "${EXPECTED_DIAGNOSTIC}" diagnostic_position)
if(diagnostic_position EQUAL -1)
  message(FATAL_ERROR "Expected diagnostic not found: ${EXPECTED_DIAGNOSTIC}\n${output}")
endif()

message(STATUS "Observed expected failure: ${EXPECTED_DIAGNOSTIC}")
