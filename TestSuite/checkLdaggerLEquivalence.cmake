if(NOT DEFINED CHECK_VALUE_EXECUTABLE)
  message(FATAL_ERROR "CHECK_VALUE_EXECUTABLE is required")
endif()

if(NOT DEFINED KRON_OUTPUT_FILE)
  message(FATAL_ERROR "KRON_OUTPUT_FILE is required")
endif()

if(NOT DEFINED STORED_OUTPUT_FILE)
  message(FATAL_ERROR "STORED_OUTPUT_FILE is required")
endif()

if(NOT DEFINED TOLERANCE)
  message(FATAL_ERROR "TOLERANCE is required")
endif()

function(read_final_lowest_eigenvalue filename result)
  if(NOT EXISTS "${filename}")
    message(FATAL_ERROR "Output file does not exist: ${filename}")
  endif()

  file(READ "${filename}" output)
  string(REGEX MATCHALL "Found lowest eigenvalue= *[-+0-9.eE]+" matches "${output}")
  if(NOT matches)
    message(FATAL_ERROR "No lowest eigenvalue found in ${filename}")
  endif()

  list(GET matches -1 last_match)
  string(REGEX REPLACE "^.*= *" "" value "${last_match}")
  set(${result}
      "${value}"
      PARENT_SCOPE)
endfunction()

read_final_lowest_eigenvalue("${KRON_OUTPUT_FILE}" kron_value)
read_final_lowest_eigenvalue("${STORED_OUTPUT_FILE}" stored_value)

execute_process(
  COMMAND "${CHECK_VALUE_EXECUTABLE}" "Found lowest eigenvalue=\\s+([.0-9eE+-]+)" "${stored_value}" "${TOLERANCE}"
          "${KRON_OUTPUT_FILE}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE stderr)
if(NOT result EQUAL 0)
  message(
    FATAL_ERROR
      "LdaggerL results differ: Kron=${kron_value}, Stored=${stored_value}, tolerance=${TOLERANCE}\n${stdout}${stderr}")
endif()

message(STATUS "LdaggerL results agree: Kron=${kron_value}, Stored=${stored_value}\n${stdout}")
