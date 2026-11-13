if(NOT DEFINED DMRG_EXECUTABLE)
  message(FATAL_ERROR "DMRG_EXECUTABLE is required")
endif()

if(NOT DEFINED INPUT211_FILE)
  message(FATAL_ERROR "INPUT211_FILE is required")
endif()

if(NOT DEFINED INPUT212_FILE)
  message(FATAL_ERROR "INPUT212_FILE is required")
endif()

if(NOT DEFINED WORK_ROOT)
  message(FATAL_ERROR "WORK_ROOT is required")
endif()

file(REMOVE_RECURSE "${WORK_ROOT}")
file(MAKE_DIRECTORY "${WORK_ROOT}/seed" "${WORK_ROOT}/resume" "${WORK_ROOT}/terminal")

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${INPUT211_FILE}"
  WORKING_DIRECTORY "${WORK_ROOT}/seed"
  RESULT_VARIABLE input211_result
  OUTPUT_VARIABLE input211_stdout
  ERROR_VARIABLE input211_stderr)
if(NOT input211_result EQUAL 0)
  message(FATAL_ERROR "Explicit-restart seed failed:\n${input211_stdout}${input211_stderr}")
endif()

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${INPUT212_FILE}"
  WORKING_DIRECTORY "${WORK_ROOT}/seed"
  RESULT_VARIABLE bootstrap_result
  OUTPUT_VARIABLE bootstrap_stdout
  ERROR_VARIABLE bootstrap_stderr)
if(NOT bootstrap_result EQUAL 0)
  message(FATAL_ERROR "Explicit-restart bootstrap failed:\n${bootstrap_stdout}${bootstrap_stderr}")
endif()

set(seed_log "${WORK_ROOT}/seed/runForinput212.cout")
file(READ "${seed_log}" seed_log_output)
if(NOT seed_log_output MATCHES "Read from file data211.hd5 succeeded")
  message(FATAL_ERROR "Bootstrap did not read data211.hd5")
endif()
if(NOT seed_log_output MATCHES "Finite loop number 0 with l=6")
  message(FATAL_ERROR "Bootstrap did not begin at finite loop zero")
endif()
if(seed_log_output MATCHES "RecoverySave: recovering from")
  message(FATAL_ERROR "Bootstrap incorrectly selected a RecoverySave file")
endif()

set(nonterminal_recovery "${WORK_ROOT}/seed/Recovery23data212.hd5")
set(terminal_recovery "${WORK_ROOT}/seed/Recovery42data212.hd5")
foreach(recovery_file IN ITEMS "${nonterminal_recovery}" "${terminal_recovery}")
  if(NOT EXISTS "${recovery_file}")
    message(FATAL_ERROR "Expected retained recovery file is missing: ${recovery_file}")
  endif()
endforeach()

file(COPY "${WORK_ROOT}/seed/data211.hd5" "${nonterminal_recovery}" DESTINATION "${WORK_ROOT}/resume")
execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${INPUT212_FILE}"
  WORKING_DIRECTORY "${WORK_ROOT}/resume"
  RESULT_VARIABLE resume_result
  OUTPUT_VARIABLE resume_stdout
  ERROR_VARIABLE resume_stderr)
set(resume_output "${resume_stdout}${resume_stderr}")
if(NOT resume_result EQUAL 0)
  message(FATAL_ERROR "RecoverySave continuation failed:\n${resume_output}")
endif()

set(resume_log "${WORK_ROOT}/resume/runForinput212.cout")
file(READ "${resume_log}" resume_log_output)
if(NOT resume_log_output MATCHES "RecoverySave: recovering from Recovery23data212.hd5")
  message(FATAL_ERROR "Continuation did not report Recovery23data212.hd5")
endif()
if(NOT resume_log_output MATCHES "Finite loop number 24 with l=6")
  message(FATAL_ERROR "Continuation did not begin at the saved loop index")
endif()
if(resume_log_output MATCHES "Finite loop number 0 with l=6")
  message(FATAL_ERROR "Continuation restarted from the bootstrap checkpoint")
endif()

file(COPY "${WORK_ROOT}/seed/data211.hd5" "${terminal_recovery}" DESTINATION "${WORK_ROOT}/terminal")
execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${INPUT212_FILE}"
  WORKING_DIRECTORY "${WORK_ROOT}/terminal"
  RESULT_VARIABLE terminal_result
  OUTPUT_VARIABLE terminal_stdout
  ERROR_VARIABLE terminal_stderr)
set(terminal_output "${terminal_stdout}${terminal_stderr}")
if(NOT terminal_result EQUAL 0)
  message(FATAL_ERROR "Terminal RecoverySave continuation failed:\n${terminal_output}")
endif()
if(NOT terminal_output MATCHES "RecoverySave: checkpoint already completed; exiting\\.")
  message(FATAL_ERROR "Terminal RecoverySave continuation did not perform a no-op:\n${terminal_output}")
endif()

message(STATUS "Explicit restart bootstraps and RecoverySave continuation resumes its output")
