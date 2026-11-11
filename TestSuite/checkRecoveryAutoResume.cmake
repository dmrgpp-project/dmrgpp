if(NOT DEFINED DMRG_EXECUTABLE)
  message(FATAL_ERROR "DMRG_EXECUTABLE is required")
endif()

if(NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE is required")
endif()

if(NOT DEFINED WORK_ROOT)
  message(FATAL_ERROR "WORK_ROOT is required")
endif()

if(NOT DEFINED CHECK_VALUE_EXECUTABLE)
  message(FATAL_ERROR "CHECK_VALUE_EXECUTABLE is required")
endif()

file(REMOVE_RECURSE "${WORK_ROOT}")
file(MAKE_DIRECTORY "${WORK_ROOT}/seed" "${WORK_ROOT}/resume" "${WORK_ROOT}/explicit" "${WORK_ROOT}/clean")

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${INPUT_FILE}"
  WORKING_DIRECTORY "${WORK_ROOT}/seed"
  RESULT_VARIABLE seed_result
  OUTPUT_VARIABLE seed_stdout
  ERROR_VARIABLE seed_stderr)
if(NOT seed_result EQUAL 0)
  message(FATAL_ERROR "RecoverySave seed run failed:\n${seed_stdout}${seed_stderr}")
endif()

set(recovery_file "${WORK_ROOT}/seed/Recovery23data210.hd5")
if(NOT EXISTS "${recovery_file}")
  message(FATAL_ERROR "Expected retained nonterminal recovery file is missing: ${recovery_file}")
endif()
file(COPY "${recovery_file}" DESTINATION "${WORK_ROOT}/resume")
file(COPY "${recovery_file}" DESTINATION "${WORK_ROOT}/explicit")
file(COPY "${WORK_ROOT}/seed/Recovery24data210.hd5" DESTINATION "${WORK_ROOT}/explicit")

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${INPUT_FILE}"
  WORKING_DIRECTORY "${WORK_ROOT}/resume"
  RESULT_VARIABLE resume_result
  OUTPUT_VARIABLE resume_stdout
  ERROR_VARIABLE resume_stderr)
set(resume_output "${resume_stdout}${resume_stderr}")
if(NOT resume_result EQUAL 0)
  message(FATAL_ERROR "Automatic RecoverySave resume failed:\n${resume_output}")
endif()
if(resume_output MATCHES "autoNumberAfterFirst: Internal error")
  message(FATAL_ERROR "Automatic RecoverySave resume hit the @auto invariant:\n${resume_output}")
endif()

file(READ "${INPUT_FILE}" explicit_input)
string(REPLACE "SolverOptions=twositedmrg,recoveryenableread;" "SolverOptions=twositedmrg,recoveryenableread,restart;"
               explicit_input "${explicit_input}")
string(APPEND explicit_input "\nRestartFilename=\"Recovery23data210\";\n")
set(explicit_input_file "${WORK_ROOT}/explicit/input210-explicit-restart.ain")
file(WRITE "${explicit_input_file}" "${explicit_input}")
execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${explicit_input_file}"
  WORKING_DIRECTORY "${WORK_ROOT}/explicit"
  RESULT_VARIABLE explicit_result
  OUTPUT_VARIABLE explicit_stdout
  ERROR_VARIABLE explicit_stderr)
if(NOT explicit_result EQUAL 0)
  message(FATAL_ERROR "Explicit restart control failed:\n${explicit_stdout}${explicit_stderr}")
endif()

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${INPUT_FILE}"
  WORKING_DIRECTORY "${WORK_ROOT}/clean"
  RESULT_VARIABLE clean_result
  OUTPUT_VARIABLE clean_stdout
  ERROR_VARIABLE clean_stderr)
if(NOT clean_result EQUAL 0)
  message(FATAL_ERROR "Clean reference run failed:\n${clean_stdout}${clean_stderr}")
endif()

set(resume_log "${WORK_ROOT}/resume/runForinput210.cout")
set(explicit_log "${WORK_ROOT}/explicit/runForinput210-explicit-restart.cout")
set(clean_log "${WORK_ROOT}/clean/runForinput210.cout")
foreach(log_file IN ITEMS "${resume_log}" "${explicit_log}" "${clean_log}")
  if(NOT EXISTS "${log_file}")
    message(FATAL_ERROR "Expected solver log is missing: ${log_file}")
  endif()
endforeach()

file(READ "${explicit_log}" explicit_output)
if(NOT explicit_output MATCHES "Finite loop number 24 with l=6")
  message(FATAL_ERROR "Explicit restart did not use RestartFilename=Recovery23data210")
endif()

file(READ "${resume_log}" resume_log_output)
if(NOT resume_log_output MATCHES "RecoverySave: recovering from Recovery23data210.hd5")
  message(FATAL_ERROR "Automatic RecoverySave resume did not report its recovery file in ${resume_log}")
endif()

file(READ "${clean_log}" clean_output)
string(REGEX MATCHALL "Found lowest eigenvalue= *[-+0-9.eE]+" clean_matches "${clean_output}")
if(NOT clean_matches)
  message(FATAL_ERROR "No lowest eigenvalue found in ${clean_log}")
endif()
list(GET clean_matches -1 clean_match)
string(REGEX REPLACE "^.*= *" "" clean_value "${clean_match}")

execute_process(
  COMMAND "${CHECK_VALUE_EXECUTABLE}" "Found lowest eigenvalue=\\s+([.0-9eE+-]+)" "${clean_value}" "2e-5"
          "${resume_log}"
  RESULT_VARIABLE compare_result
  OUTPUT_VARIABLE compare_stdout
  ERROR_VARIABLE compare_stderr)
if(NOT compare_result EQUAL 0)
  message(FATAL_ERROR "Recovered and clean energies differ:\n${compare_stdout}${compare_stderr}")
endif()

execute_process(
  COMMAND "${CHECK_VALUE_EXECUTABLE}" "Found lowest eigenvalue=\\s+([.0-9eE+-]+)" "${clean_value}" "2e-5"
          "${explicit_log}"
  RESULT_VARIABLE explicit_compare_result
  OUTPUT_VARIABLE explicit_compare_stdout
  ERROR_VARIABLE explicit_compare_stderr)
if(NOT explicit_compare_result EQUAL 0)
  message(
    FATAL_ERROR "Explicit restart and clean energies differ:\n${explicit_compare_stdout}${explicit_compare_stderr}")
endif()

message(STATUS "Automatic RecoverySave resume matches the clean run")
