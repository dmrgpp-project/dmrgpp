# Source this file from a fresh Frontier shell; it changes the current module stack.
module reset
module load PrgEnv-cray
module load craype/2.7.36
module swap cce cce/21.0.0
module load cray-mpich/9.1.0
module load rocm/6.4.2
module load craype-accel-amd-gfx90a
module load cmake
module load gsl
module load darshan-runtime/3.4.7-mpi

module is-loaded darshan-runtime/3.4.7-mpi || {
  print -u2 "Darshan runtime did not activate"
  return 1
}

print -r -- "CCE: $(command -v CC)"
CC --version | head -n 1
print -r -- "MPICH: $CRAY_MPICH_DIR"
module -t list
