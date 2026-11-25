find_program(PDFLATEX_EXE pdflatex REQUIRED)
find_program(BIBTEX_EXE bibtex REQUIRED)
find_program(PERL_EXE perl REQUIRED)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/manual.tex
  COMMAND bash -c "cp ${CMAKE_CURRENT_SOURCE_DIR}/*.pdf ${CMAKE_CURRENT_BINARY_DIR}/"
  COMMAND bash -c "cp ${CMAKE_CURRENT_SOURCE_DIR}/*.png ${CMAKE_CURRENT_BINARY_DIR}/"
  COMMAND
    find ${dmrgpp_SOURCE_DIR}/dmrg -iname "*.h" -or -iname "*.cpp" | ${PERL_EXE}
    ${dmrgpp_SOURCE_DIR}/PsimagLite/scripts/doc.pl ${CMAKE_CURRENT_SOURCE_DIR}/manual.ptex ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating manual.tex from manual.ptex"
  #  VERBATIM
)

add_custom_target(
  manual.pdf ALL
  COMMAND ${PDFLATEX_EXE} manual.tex
  COMMAND ${BIBTEX_EXE} manual.aux
  COMMAND ${PDFLATEX_EXE} manual.tex
  COMMAND ${PDFLATEX_EXE} manual.tex
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating manual.pdf"
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/manual.tex
  VERBATIM)
