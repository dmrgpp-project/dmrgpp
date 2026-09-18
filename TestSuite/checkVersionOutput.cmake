if(NOT DEFINED DMRG_EXECUTABLE)
  message(FATAL_ERROR "DMRG_EXECUTABLE is required")
endif()

if(NOT DEFINED EXPECTED_VERSION)
  message(FATAL_ERROR "EXPECTED_VERSION is required")
endif()

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -V
  RESULT_VARIABLE result
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE stderr)

set(output "${stdout}${stderr}")
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Version command failed:\n${output}")
endif()

string(ASCII 27 escape)
set(expected "DMRG++::dmrg${escape}[38;5;240m [master ${EXPECTED_VERSION}] ${escape}[0m")
string(FIND "${output}" "${expected}" expected_position)
if(expected_position EQUAL -1)
  message(FATAL_ERROR "Expected version logo not found:\n${expected}\nActual output:\n${output}")
endif()

message(STATUS "Observed expected version logo")
