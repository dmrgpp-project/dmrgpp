#!/usr/bin/env bash
#SBATCH --job-name=cincuenta-star
#SBATCH --output=cincuenta-star-%j.out
#SBATCH --error=cincuenta-star-%j.err

set -euo pipefail

. /lustre/orion/cph102/proj-shared/epd/dmrgpp/build-aux/frontier_load_modules.sh

# Override these with sbatch --export or edit a private cluster-specific copy.
: "${CASE:=nb7_u2p5}"
: "${RANKS:=${SLURM_NTASKS:-16}}"


export CINCUENTA_EXE=/lustre/orion/cph102/proj-shared/epd/dmrgpp/build_cce21_mpich91_hip/cincuenta/src/cincuenta
export SRUN_ARGS='--ntasks-per-node=8 --cpus-per-task=7 --threads-per-core=1 --distribution=block:block --cpu-bind=cores --gpus-per-task=1'
export RUN_ROOT=$(pwd)
script_dir=$(pwd)/../cincuenta/jobs/garcia2004_star
exec "$script_dir/run_case.sh" "$CASE" "$RANKS"
