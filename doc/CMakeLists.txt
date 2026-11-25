find_program(PDFLATEX_EXE pdflatex REQUIRED)
find_program(PERL_EXE perl REQUIRED)

file(
  GLOB_RECURSE
  MANUAL_SOURCE_FILES
  CONFIGURE_DEPENDS
  ${dmrgpp_SOURCE_DIR}/dmrg/*.h
  ${dmrgpp_SOURCE_DIR}/dmrg/*.cpp
  ${dmrgpp_SOURCE_DIR}/LanczosPlusPlus/*.h
  ${dmrgpp_SOURCE_DIR}/LanczosPlusPlus/*.cpp
  ${dmrgpp_SOURCE_DIR}/PsimagLite/*.h
  ${dmrgpp_SOURCE_DIR}/PsimagLite/*.cpp
  ${dmrgpp_SOURCE_DIR}/util/*.h
  ${dmrgpp_SOURCE_DIR}/util/*.cpp)
list(JOIN MANUAL_SOURCE_FILES "\n" MANUAL_SOURCE_LIST)
file(
  GENERATE
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/manual-sources.txt
  CONTENT "${MANUAL_SOURCE_LIST}\n")

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/manual.tex
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/wft-eps-converted-to.pdf
          ${CMAKE_CURRENT_BINARY_DIR}/wft.pdf
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/dmrgV2LogoBW.png
          ${CMAKE_CURRENT_BINARY_DIR}/dmrgV2LogoBW.png
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/codemap.png
          ${CMAKE_CURRENT_BINARY_DIR}/codemap.png
  COMMAND
    ${CMAKE_COMMAND} -DPERL_EXE=${PERL_EXE} -DDOC_SCRIPT=${dmrgpp_SOURCE_DIR}/PsimagLite/scripts/doc.pl
    -DMANUAL_INPUT=${CMAKE_CURRENT_SOURCE_DIR}/manual.ptex -DDOC_ROOT=${CMAKE_CURRENT_SOURCE_DIR}
    -DOUTPUT_DIR=${CMAKE_CURRENT_BINARY_DIR} -DSOURCE_LIST=${CMAKE_CURRENT_BINARY_DIR}/manual-sources.txt -P
    ${CMAKE_CURRENT_SOURCE_DIR}/GenerateManual.cmake
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating manual.tex from manual.ptex"
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/GenerateManual.cmake
          ${CMAKE_CURRENT_SOURCE_DIR}/manual.ptex
          ${CMAKE_CURRENT_SOURCE_DIR}/wft-eps-converted-to.pdf
          ${CMAKE_CURRENT_SOURCE_DIR}/dmrgV2LogoBW.png
          ${CMAKE_CURRENT_SOURCE_DIR}/codemap.png
          ${dmrgpp_SOURCE_DIR}/LICENSE
          ${dmrgpp_SOURCE_DIR}/NOTICE
          ${dmrgpp_SOURCE_DIR}/THIRD_PARTY_NOTICES
          ${dmrgpp_SOURCE_DIR}/PsimagLite/scripts/doc.pl
          ${CMAKE_CURRENT_BINARY_DIR}/manual-sources.txt
          ${MANUAL_SOURCE_FILES}
  VERBATIM)

set(MANUAL_BIBTEX_COMMAND)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/thesis.bib)
  find_program(BIBTEX_EXE bibtex REQUIRED)
  list(
    APPEND
    MANUAL_BIBTEX_COMMAND
    COMMAND
    ${CMAKE_COMMAND}
    -E
    copy_if_different
    ${CMAKE_CURRENT_SOURCE_DIR}/thesis.bib
    ${CMAKE_CURRENT_BINARY_DIR}/thesis.bib
    COMMAND
    ${BIBTEX_EXE}
    manual.aux)
endif()

add_custom_target(
  manual.pdf ALL
  COMMAND ${PDFLATEX_EXE} manual.tex ${MANUAL_BIBTEX_COMMAND}
  COMMAND ${PDFLATEX_EXE} manual.tex
  COMMAND ${PDFLATEX_EXE} manual.tex
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating manual.pdf"
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/manual.tex
  VERBATIM)
