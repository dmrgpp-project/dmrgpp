if(NOT DEFINED DMRG_EXECUTABLE
   OR NOT DEFINED INPUT_FILE
   OR NOT DEFINED WORKING_DIRECTORY)
  message(FATAL_ERROR "DMRG_EXECUTABLE, INPUT_FILE, and WORKING_DIRECTORY are required")
endif()

file(REMOVE "${WORKING_DIRECTORY}/data31_operator.hd5" "${WORKING_DIRECTORY}/data31_operator00.hd5"
     "${WORKING_DIRECTORY}/kron_data31_operator_0.txt" "${WORKING_DIRECTORY}/runFordata31_operator.cout")

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${INPUT_FILE}"
  WORKING_DIRECTORY "${WORKING_DIRECTORY}"
  RESULT_VARIABLE dmrg_result
  OUTPUT_FILE "${WORKING_DIRECTORY}/input31_kron_operator.stdout"
  ERROR_FILE "${WORKING_DIRECTORY}/input31_kron_operator.stderr")

if(NOT dmrg_result EQUAL 0)
  message(FATAL_ERROR "dmrg failed with exit code ${dmrg_result}")
endif()

if(NOT EXISTS "${WORKING_DIRECTORY}/kron_data31_operator_0.txt")
  message(FATAL_ERROR "dmrg did not produce kron_data31_operator_0.txt")
endif()
