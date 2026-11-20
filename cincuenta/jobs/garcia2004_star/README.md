# Garcia--Hallberg--Rozenberg 2004 star-DMFT jobs

These are cluster workloads for comparing cincuenta's finite **star bath** with
the equilibrium Bethe-DMFT results of:

> D. J. Garcia, K. Hallberg, and M. J. Rozenberg, *Dynamical Mean Field Theory
> with the Density Matrix Renormalization Group*, Phys. Rev. Lett. 93, 246403
> (2004), [arXiv:cond-mat/0403169](https://arxiv.org/abs/cond-mat/0403169v1).

They are physics-motivated scale tests, not an exact reproduction. The paper
uses two continued-fraction chains with as many as 45 total sites. Cincuenta
fits a discrete star bath on the Matsubara axis. Equal chain and star site
counts therefore are not equivalent finite-size approximations; compare
convergence of Green functions, low-energy weight, and gaps instead.

## Cases and suggested order

| Case | Purpose | Bath/total sites | Kept states | Matsubaras | Real grid | Iterations | Default ranks |
|---|---|---:|---:|---:|---:|---:|---:|
| `nb7_u1p5` | low-interaction comparison | 7/8 | 128 | 32 positive | 129 | 8 | 16 |
| `nb7_u2p5` | pilot | 7/8 | 128 | 32 positive | 129 | 8 | 16 |
| `nb11_u2p0` | metallic control | 11/12 | 300 | 64 positive | 257 | 15 | 32 |
| `nb11_u3p5` | insulating control | 11/12 | 300 | 64 positive | 257 | 15 | 32 |
| `nb11_u2p5` | primary near-transition run | 11/12 | 300 | 64 positive | 257 | 15 | 32 |
| `nb15_u2p5` | large stress run | 15/16 | 600 | 128 positive | 513 | 20 | 64 |

Run the pilot first. Then run the two controls before interpreting the primary
`U=2.5` result. Attempt the 15-bath case only after memory and elapsed time are
known for 11 bath sites.

All inputs use:

- Bethe half-bandwidth `D=1`, represented by
  `LatticeGf="energy,semicircular,2"`;
- `FicticiousBeta=64`;
- the explicitly centered interaction convention, `ChemicalPotential=0`,
  consistent with the paper and with `ModelParams` placing `-U/2` on the
  impurity;
- equal up/down particle sectors at half filling;
- `FitOptions=particleholesymmetric`, which enforces mirrored bath energies
  and hybridizations around a zero-energy bath site;
- a constrained initial star bath whose expanded hybridizations satisfy
  `sum(V_i^2)=1/4=t^2`;
- five correction-vector finite sweeps. The former three-sweep schedule gave
  nonreproducible Nb7 DMFT updates; this setting was validated by repeated
  8- and 16-rank Frontier GPU runs.

These jobs require the centered-equilibrium particle-hole fix merged into this
branch. In that corrected convention, the impurity Hamiltonian already contains
the `-U/2` shift and the fit therefore requires `ChemicalPotential=0`. Each
`InitBathVector` contains `(Nb+1)/2` independent hybridizations followed by
`(Nb-1)/2` positive bath energies; cincuenta expands it to `2*Nb` scalar
parameters describing `Nb` symmetric bath sites.

Set
`CINCUENTA_EXE=/absolute/path/to/cincuenta` if the cluster build is elsewhere.
Do not invoke an MPI-linked executable directly; even a one-rank run should go
through `srun` or `mpiexec`.

## Interactive or existing allocation

`run_case.sh` detects Slurm and uses `srun`; elsewhere it uses `mpiexec`.
It relies on the site's configured `srun` MPI integration; set the site-required
`SLURM_MPI_TYPE` (for example, `pmix`) in the job environment when necessary.
Every invocation creates a timestamped directory under
`runs/garcia2004_star/` and records the commit, executable/input checksums,
launcher, environment, elapsed time, and exit status.

```bash
# In a 16-rank Slurm allocation:
cincuenta/jobs/garcia2004_star/run_case.sh nb7_u2p5 16

# In a 32-rank allocation, with explicit Frontier CPU placement:
SRUN_ARGS='--distribution=block:block --cpu-bind=cores' \
  cincuenta/jobs/garcia2004_star/run_case.sh nb11_u2p5 32

# One rank per Frontier GPU/GCD on a single node:
SRUN_ARGS='--ntasks-per-node=8 --cpus-per-task=1 --distribution=block:block --cpu-bind=cores --gpus-per-task=1 --gpu-bind=closest' \
  cincuenta/jobs/garcia2004_star/run_case.sh nb7_u2p5 8

# Non-Slurm MPI, with an out-of-tree executable:
CINCUENTA_EXE=/path/to/build/cincuenta/src/cincuenta \
  cincuenta/jobs/garcia2004_star/run_case.sh nb11_u2p0 32
```

Useful overrides are documented by:

```bash
cincuenta/jobs/garcia2004_star/run_case.sh --help
```

In particular, `RUN_ROOT` can point at a large scratch filesystem and
`SRUN_ARGS` can supply site-specific CPU/GPU placement options. `SRUN_ARGS` is
split on whitespace, so its arguments cannot themselves contain spaces. The
fully expanded launch command is recorded in `manifest.txt`. Keep each run in
its own directory: DMRG restart and frequency filenames are not safe to share
between simultaneous cases.

### Local CPU comparison: `U=1.5`, `U=2.5`, and `U=3.5`

Use the dedicated `nb7_u1p5_cpu`, `nb7_u2p5_cpu`, and `nb7_u3p5_cpu` inputs
for CPU builds. These differ only in `HubbardU` and `RootOutputname`.
They retain the original GPU-oriented inputs unchanged and remove only
`BatchedGemm` from `SolverOptions` (the CPU BatchedGemm path requires dense
matrices); this execution-only setting selects the standard matvec and does
not change the physics. Run each case with 8 MPI ranks and one thread per rank:

```bash
export CINCUENTA_EXE=/absolute/path/to/build_mpi/cincuenta/src/cincuenta
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cincuenta/jobs/garcia2004_star/run_case.sh nb7_u1p5_cpu
cincuenta/jobs/garcia2004_star/run_case.sh nb7_u2p5_cpu
cincuenta/jobs/garcia2004_star/run_case.sh nb7_u3p5_cpu
```

`U=3.5` supplies a high-interaction control; this three-point comparison is
not a determination of the critical interaction. See
[`CPU_RESULTS.md`](CPU_RESULTS.md) for the local diagnostic results and
execution-path caveats.

### Exact diagonalization comparison

The `nb7_u1p5_ed`, `nb7_u2p5_ed`, and `nb7_u3p5_ed` inputs match their CPU DMRG
counterparts except for the solver and output name. Run with the MPI-enabled
`build_mpi` executable; ED defaults to one rank because it has no distributed
frequency workload:

```bash
export CINCUENTA_EXE=/absolute/path/to/build_mpi/cincuenta/src/cincuenta
cincuenta/jobs/garcia2004_star/run_case.sh nb7_u1p5_ed
cincuenta/jobs/garcia2004_star/run_case.sh nb7_u2p5_ed
cincuenta/jobs/garcia2004_star/run_case.sh nb7_u3p5_ed
```

These matched comparisons retain iteration cap 8, `FicticiousBeta=64`, and
`Matsubaras=32` (as well as the DMRG-specific input keys) for direct provenance.
**The original repeats exposed ED spectral-weight and normalization defects.
The shared-engine factor-four defect is now fixed and the cases have been
rerun, but the finite-tail normalization policy remains unchanged. See
[`ED_RESULTS.md`](ED_RESULTS.md) before interpreting ED amplitudes or comparing
these solutions quantitatively with DMRG.**

## Slurm batch submission

`slurm_job.sh` contains deliberately minimal portable directives. Supply the
site-specific account, partition, node count, memory, and wall time on the
`sbatch` command line or in a private copy. Command-line options override the
16-task defaults in the file.

```bash
# Pilot
sbatch -A cph102 -N2 --ntasks=16 --time=00:30:00 \
--export=ALL,CASE=nb7_u2p5,RANKS=16 -p service -q develop ./slurm_job.sh

# Primary run
sbatch -A cph102 -N4 --ntasks=32 --time=00:30:00 \
--export=ALL,CASE=nb11_u2p5,RANKS=32 -p batch -q debug ./slurm_job.sh

```

For multiple nodes, request enough nodes for the task count according to the
site's tasks-per-node policy. The script deliberately does not guess this.

## Outputs and quick comparison

Important files in each run directory are:

- `manifest.txt`: provenance, launch command, elapsed time, and status;
- `cincuenta.log`: rank-zero solver output, final bath parameters, and DMFT
  iteration errors;
- `gimp_<solver>.txt`: final Matsubara impurity Green function;
- `latticeG_<solver>.txt`: final Matsubara local lattice Green function;
- `gimp_<solver>_real.txt`: broadened real-axis impurity Green function;

Here `<solver>` is read from `ImpuritySolver` in the run input (`dmrg` is
assumed for older inputs without that key).
- `launcher.log`: complete, unfiltered launcher and stderr output.

The interactive launcher suppresses only known high-volume benign DMRG
restart diagnostics, while retaining them in `launcher.log`. Rank-zero
`DMFT_PROGRESS` lines identify each SCF iteration and the final real-axis
solve. Rank-zero `DMFT_TIMING` lines report the GS-run wall time, distributed
correction-vector runner time (slowest rank and summed work), and rank-zero
text-output processing time for each Matsubara or real-axis solve. Set
`SHOW_ALL_LAUNCHER_OUTPUT=1` to show the unfiltered stream in the terminal
when diagnosing a failure.

Summarize one run or every child of the run root:

```bash
python3 cincuenta/jobs/garcia2004_star/summarize_runs.py \
  runs/garcia2004_star
```

The summary marks a run complete only when its manifest reports a successful
exit and all three Green-function tables have their expected row counts. It
returns a nonzero status if any selected run is incomplete; use
`--allow-incomplete` only when monitoring work still in progress.

The reported `A_nearest_zero` is `-Im G(omega)/pi` at the real-grid point
nearest zero. Treat it as a broadened finite-bath diagnostic, not a direct
critical-order parameter. The more robust initial comparison is the low-
frequency shape of `-Im G(i omega_n)` between `U=2.0`, `2.5`, and `3.5`.

Plot both Green-function axes (the ordinate is `-Im G`, not divided by π):

```bash
python3 cincuenta/jobs/garcia2004_star/plot_green_functions.py \
  build_mpi/runs/garcia2004_star_cpu/nb7_* \
  --output build_mpi/plots/nb7-green-functions.svg
```

## Scaling interpretation and stopping criteria

MPI distributes independent correction-vector frequencies. Rank 0 still runs
each ground-state DMRG and assembles the result, so frequency scaling cannot
remove that serial cost. Using more ranks than useful independent frequencies
will not help.

Before advancing to the next case, require:

1. successful exit and finite three-column Green-function files;
2. decreasing DMFT errors, preferably convergence below `1e-4`;
3. approximately particle-hole-symmetric final bath parameters and spectra;
4. physically distinct metallic (`U=2.0`) and insulating (`U=3.5`) controls;
5. acceptable truncation behavior and no memory pressure in the DMRG log.

If a case reaches its iteration limit without convergence, retain the full run
directory. Do not silently increase bath size or interpret the real-axis line
shape before inspecting the fit and truncation diagnostics.
