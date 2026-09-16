# SPDX-FileCopyrightText: Copyright (c) 2026, UT-Battelle, LLC
# SPDX-License-Identifier: GPL-3.0-or-later

foreach(required_variable IN ITEMS PERL_EXE DOC_SCRIPT MANUAL_INPUT DOC_ROOT OUTPUT_DIR SOURCE_LIST)
  if(NOT DEFINED ${required_variable})
    message(FATAL_ERROR "GenerateManual.cmake requires ${required_variable}")
  endif()
endforeach()

execute_process(COMMAND "${PERL_EXE}" "${DOC_SCRIPT}" "${MANUAL_INPUT}" "${DOC_ROOT}" "${OUTPUT_DIR}"
                INPUT_FILE "${SOURCE_LIST}" COMMAND_ERROR_IS_FATAL ANY)

file(READ "${OUTPUT_DIR}/manual.tex" generated_manual)
if(generated_manual MATCHES "ERROR: Label not found")
  message(FATAL_ERROR "Generated manual.tex contains unresolved PsimagDoc labels")
endif()
