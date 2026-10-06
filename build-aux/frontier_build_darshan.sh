# Source frontier_load_modules_darshan.sh first, then run this from the source tree.
# The separate Spack environment supplies HDF5 and Boost only. Deactivate it
# before running this script so CMake sees the loaded Cray wrappers unchanged.

spack_env=${FRONTIER_DARSHAN_SPACK_ENV:-/ccs/home/pdoak/spack_env/dmrgpp-cce2100-darshan}
repo="$(git rev-parse --show-toplevel)"
build=${FRONTIER_DARSHAN_BUILD_DIR:-"$repo/build_cce2100_mpich910_hip_darshan"}

export HDF5_ROOT="$(spack -e "$spack_env" location -i hdf5)"
export Boost_ROOT="$(spack -e "$spack_env" location -i boost)"

print -r -- "CXX wrapper: $(command -v CC)"
CC --version | head -n 1
print -r -- "MPICH: $CRAY_MPICH_DIR"
print -r -- "HDF5: $HDF5_ROOT"
print -r -- "Boost: $Boost_ROOT"

[[ "$CRAY_MPICH_DIR" == */mpich/9.1.0/* ]] || {
  print -u2 "Unexpected Cray MPICH flavor: $CRAY_MPICH_DIR"
  return 1
}
[[ -f "$CRAY_MPICH_DIR/include/mpi.h" ]] || {
  print -u2 "Missing MPI header"
  return 1
}
[[ -f "$HDF5_ROOT/include/H5Cpp.h" ]] || {
  print -u2 "Missing HDF5 C++ header"
  return 1
}
[[ -f "$Boost_ROOT/include/boost/version.hpp" ]] || {
  print -u2 "Missing Boost headers"
  return 1
}
[[ ! -e "$build/CMakeCache.txt" ]] || {
  print -u2 "Refusing to reuse existing build: $build"
  return 1
}

module -t list 2> "$repo/frontier-cce2100-darshan.modules"

cmake -S "$repo" -B "$build" \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_C_COMPILER="$(command -v cc)" \
  -DCMAKE_CXX_COMPILER="$(command -v CC)" \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DBUILD_TESTING=ON \
  -DDMRGPP_USE_MPI=ON \
  -DDMRG_BATCHED_GEMM_TYPE=Kokkos \
  -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ARCH_AMD_GFX90A=ON \
  -DBLA_VENDOR=All \
  -DHDF5_ROOT="$HDF5_ROOT" \
  -DBoost_ROOT="$Boost_ROOT" \
  -Dhip_DIR="$ROCM_PATH/lib/cmake/hip" \
  -DCMAKE_PREFIX_PATH="$HDF5_ROOT;$Boost_ROOT;$ROCM_PATH"
