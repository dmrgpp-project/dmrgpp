#!/usr/bin/env bash
# Create the isolated CCE 21.0.0 / Cray MPICH 9.1.0 dependency environment.
# Source frontier_load_modules_darshan.sh first. This script neither installs
# packages nor modifies the known-good dmrgpp Spack environment.
set -euo pipefail

spack_env=${FRONTIER_DARSHAN_SPACK_ENV:-/ccs/home/pdoak/spack_env/dmrgpp-cce2100-darshan}
install_root=${FRONTIER_DARSHAN_SPACK_INSTALL_ROOT:-/lustre/orion/cph102/proj-shared/epd/spack-installs/dmrgpp-cce2100-mpich910-darshan}

: "${CRAY_MPICH_DIR:?Load the experimental Frontier modules before creating this environment.}"
[[ $CRAY_MPICH_DIR == */mpich/9.1.0/* ]] || {
  echo "Unexpected Cray MPICH directory: $CRAY_MPICH_DIR" >&2
  exit 1
}
[[ -x /opt/cray/pe/cce/21.0.0/bin/crayCC ]] || {
  echo "CCE 21.0.0 compiler is not installed at the expected prefix" >&2
  exit 1
}
[[ ! -e $spack_env/spack.yaml ]] || {
  echo "Refusing to overwrite existing Spack environment: $spack_env" >&2
  exit 1
}

mkdir -p "$spack_env"
cat > "$spack_env/spack.yaml" <<EOF
spack:
  specs:
  - "hdf5@1.14.6 +cxx +hl +mpi ~fortran +shared +tools %cce@21.0.0 ^cray-mpich@9.1.0%cce@21.0.0"
  - "boost@1.89.0 +shared %cce@21.0.0"
  concretizer:
    unify: true
  config:
    install_tree:
      root: $install_root
      padded_length: 128
  packages:
    all:
      providers:
        mpi: [cray-mpich]
    cce:
      externals:
      - spec: cce@21.0.0
        prefix: /opt/cray/pe/cce/21.0.0
        extra_attributes:
          compilers:
            c: /opt/cray/pe/cce/21.0.0/bin/craycc
            cxx: /opt/cray/pe/cce/21.0.0/bin/crayCC
            fortran: /opt/cray/pe/cce/21.0.0/bin/crayftn
      buildable: false
    cray-mpich:
      externals:
      - spec: cray-mpich@9.1.0%cce@21.0.0
        prefix: $CRAY_MPICH_DIR
      buildable: false
EOF

printf 'Created isolated Spack environment: %s\n' "$spack_env"
printf 'Dedicated install tree: %s\n' "$install_root"
printf 'Review then concretize with:\n  spack -e %q concretize --fresh\n  spack -e %q spec -Il\n' "$spack_env" "$spack_env"
