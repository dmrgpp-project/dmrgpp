# Add -Wextra in the future
add_library(
  psimaglite
  PsimagLite/MersenneTwister.cpp
  PsimagLite/Matrix.cpp
  PsimagLite/Mpi.cpp
  PsimagLite/Concurrency.cpp
  PsimagLite/ProgressIndicator.cpp
  PsimagLite/MemResolv.cpp
  PsimagLite/PsimagLite.cpp
  PsimagLite/PsiBase64.cpp
  PsimagLite/SpecialFunctions.cpp
  PsimagLite/Io/TypeToH5.cpp
  PsimagLite/TridiagonalMatrix.cpp
  PsimagLite/PredicateSimple.cpp
  PsimagLite/Ainur/AinurSpirit.cpp
  PsimagLite/Ainur/AinurConvert.cpp)

# Resolve Kokkos and KokkosKernels as a coherent pair. An installed
# KokkosKernels package selects the Kokkos installation it was built with. If
# that pair is unavailable, build both dependencies from source.
include(FetchContent)

if(TARGET Kokkos::kokkos OR TARGET Kokkos::kokkoskernels)
  message(FATAL_ERROR "Kokkos targets must not exist before DMRG++ resolves Kokkos and KokkosKernels as a pair")
endif()

set(_dmrgpp_kokkoskernels_package_found FALSE)
if(NOT FETCHCONTENT_TRY_FIND_PACKAGE_MODE STREQUAL "NEVER"
   AND NOT FETCHCONTENT_SOURCE_DIR_KOKKOS
   AND NOT FETCHCONTENT_SOURCE_DIR_KOKKOSKERNELS)
  list(APPEND CMAKE_MESSAGE_INDENT "[KokkosKernels] ")
  find_package(KokkosKernels 5.1 QUIET CONFIG NAMES KokkosKernels)
  list(POP_BACK CMAKE_MESSAGE_INDENT)
  if(KokkosKernels_FOUND)
    set(_dmrgpp_kokkoskernels_package_found TRUE)
  endif()
endif()

if(_dmrgpp_kokkoskernels_package_found)
  # Remove values forced for a source build by an earlier configure. Newer
  # installed packages export these as normal variables that can be checked.
  unset(KokkosKernels_ENABLE_TPL_BLAS CACHE)
  unset(KokkosKernels_ENABLE_TPL_LAPACK CACHE)
  if(KokkosKernels_VERSION VERSION_GREATER_EQUAL 5.2.99)
    if(NOT KokkosKernels_ENABLE_TPL_BLAS)
      message(FATAL_ERROR "KokkosKernels must be configured with BLAS support!")
    endif()
    if(NOT KokkosKernels_ENABLE_TPL_LAPACK)
      message(FATAL_ERROR "KokkosKernels must be configured with LAPACK support!")
    endif()
  endif()
else()
  # Do not substitute either source dependency with an independently installed
  # package; the two source builds must use the same Kokkos target.
  block()
  set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE NEVER)
  FetchContent_Declare(
    Kokkos
    GIT_REPOSITORY https://github.com/kokkos/kokkos
    GIT_TAG 5.2.0
    SYSTEM)
  FetchContent_Declare(
    KokkosKernels
    GIT_REPOSITORY https://github.com/kokkos/kokkos-kernels
    GIT_TAG 5.2.0
    SYSTEM)
  endblock()

  list(APPEND CMAKE_MESSAGE_INDENT "[Kokkos] ")
  FetchContent_MakeAvailable(Kokkos)
  list(POP_BACK CMAKE_MESSAGE_INDENT)

  # Since DMRG++ requires BLAS and LAPACK, require the source build to use both
  # as TPLs.
  set(KokkosKernels_ENABLE_TPL_BLAS
      ON
      CACHE BOOL "Whether to enable the BLAS TPL" FORCE)
  set(KokkosKernels_ENABLE_TPL_LAPACK
      ON
      CACHE BOOL "Whether to enable the LAPACK TPL" FORCE)
  list(APPEND CMAKE_MESSAGE_INDENT "[KokkosKernels] ")
  FetchContent_MakeAvailable(KokkosKernels)
  list(POP_BACK CMAKE_MESSAGE_INDENT)
endif()

if(NOT TARGET Kokkos::kokkos OR NOT TARGET Kokkos::kokkoskernels)
  message(FATAL_ERROR "Kokkos and KokkosKernels did not provide their required CMake targets")
endif()

# Third-party targets should not inherit this project's clang-tidy checks.
set(KOKKOS_TARGETS kokkosalgorithms kokkoscontainers kokkoscore kokkossimd kokkoskernels)
foreach(target ${KOKKOS_TARGETS})
  if(TARGET ${target})
    set_target_properties(${target} PROPERTIES CXX_CLANG_TIDY "")
  endif()
endforeach()
if(TARGET kokkoskernels)
  set_target_properties(kokkoskernels PROPERTIES CXX_CLANG_TIDY "")
endif()

target_link_libraries(psimaglite PUBLIC Kokkos::kokkoskernels)
target_sources(psimaglite PRIVATE PsimagLite/KokkosGemm.cpp)

target_include_directories(psimaglite PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/Io
                                             ${CMAKE_CURRENT_SOURCE_DIR}/Ainur)

target_compile_features(psimaglite PUBLIC cxx_std_17)
target_link_libraries(psimaglite PUBLIC loki)

add_library(psimaglite::psimaglite ALIAS psimaglite)

if(DMRGPP_USE_MPI)
  target_compile_definitions(psimaglite PUBLIC DMRGPP_USE_MPI)
  target_link_libraries(psimaglite PUBLIC MPI::MPI_CXX)
endif()

# Prefer HDF5 parallel when MPI is enabled.
set(HDF5_PREFER_PARALLEL ${DMRGPP_USE_MPI})
find_package(
  HDF5 REQUIRED
  COMPONENTS CXX HL
  MODULE)
if(DMRGPP_USE_MPI AND NOT HDF5_IS_PARALLEL)
  message(WARNING "HDF5 that was found \(${HDF5_DIR}\) is not parallel")
endif()
target_link_libraries(psimaglite PUBLIC HDF5::HDF5)

find_package(Boost REQUIRED)
target_link_libraries(psimaglite PUBLIC Boost::boost)

find_package(Threads)
if(Threads_FOUND)
  target_compile_definitions(psimaglite PUBLIC USE_PTHREADS)
  target_link_libraries(psimaglite PUBLIC Threads::Threads)
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(psimaglite PUBLIC BLAS::BLAS)
target_link_libraries(psimaglite PUBLIC LAPACK::LAPACK)

find_package(GSL REQUIRED)
target_link_libraries(psimaglite PUBLIC GSL::gsl)

if(BUILD_SHARED_LIBS)
  install(TARGETS psimaglite LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(BUILD_TESTING)
  add_subdirectory(tests)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/PsimagLite/Ainur/function.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  add_executable(Ainur_Test1 PsimagLite/Ainur/test.cpp)
  target_link_libraries(Ainur_Test1 psimaglite::psimaglite)
  add_test(NAME Ainur_Test1 COMMAND Ainur_Test1 ${CMAKE_CURRENT_SOURCE_DIR}/PsimagLite/Ainur/prelude.ain
                                    ${CMAKE_CURRENT_SOURCE_DIR}/PsimagLite/Ainur/input0.ain)

  add_executable(Ainur_Test2 PsimagLite/Ainur/testReal.cpp)
  target_link_libraries(Ainur_Test2 psimaglite::psimaglite)
  add_test(NAME Ainur_Test2 COMMAND Ainur_Test2 ${CMAKE_CURRENT_SOURCE_DIR}/PsimagLite/Ainur/input1.ain)
endif()
