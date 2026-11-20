# Local Nb7 CPU diagnostic (2026-11-20)

## Purpose and setup

Check whether the local equilibrium calculation shows insulating behavior near
`HubbardU=1.5` with `LatticeGf="energy,semicircular,2"` (full bandwidth 2,
half-bandwidth 1). This is a small-bath diagnostic, not a critical-U estimate.

The jobs were imported from `trusted/feature/cincuenta-mpi` at
`7c87c55f305e30ed6ec140ba9fb744d72e863bd0`. Local CPU inputs use the same
Nb7 settings: 32 positive Matsubaras, fictitious beta 64, centered interaction
and zero chemical potential, particle-hole-symmetric initial bath and fit,
128 GS kept states, and five correction-vector sweeps with 256 kept states.
The DMFT iteration cap is 8 and requested error tolerance is `1e-4`.
The real-axis grid has 129 points from -4 to 4 with broadening 0.1.

Runs used 8 MPI ranks, one BLAS/OpenMP thread per rank, on the 12-CPU aarch64
host `amber_stone`. Build configuration:

```sh
cmake -S . -B build_mpi -DDMRGPP_USE_MPI=ON \
  -DDMRG_BATCHED_GEMM_TYPE=CPU -DBUILD_TESTING=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build_mpi --target cincuenta check_value -j4
```

Clang 21.1.4 and Open MPI 5.0.10 came from the active
`scientific-workstation` Spack environment. The environment view did not expose
`prterun`; `PRTE_PREFIX` was overridden per launch to the installed PRRTE
prefix, without changing the environment itself.

## CPU execution-path caveat

An initial launch of the unchanged `nb7_u1p5` input aborted during SCF iteration
0 with `FATAL: Matrix isn't dense`. `BatchedGemmCpu.h` requests dense operator
blocks, but some blocks are sparse. The failed output is preserved under
`build_mpi/runs/garcia2004_star/nb7_u1p5_20261120T181420Z`.

The successful `*_cpu` inputs remove `BatchedGemm` from `SolverOptions`, using
the standard sparse-capable matrix-vector path instead. No Hamiltonian or
convergence parameters were changed. **These runs do not validate the original
GPU BatchedGemm implementation.**

## Results

| U | SCF solves performed | Final reported error | Meets tolerance? | -Im G(i omega_0) | A(0), eta=0.1 | Elapsed seconds |
|---|---:|---:|---|---:|---:|---:|
| 1.5 | 3 | 7.53628e-6 | yes | 1.81774 | 0.531676 | 33 |
| 2.5 | 6 | 7.98930e-5 | yes | 1.42815 | 0.352261 | 56 |
| 3.5 | 8 | 4.26804e-2 | **no** | 0.247699 | 0.0495733 | 54 |

Here `omega_0 = pi/64 = 0.0490874`; `A(0) = -Im G(0+i eta)/pi`.
SCF counts above count solves, not the zero-based index printed by the solver's
final convergence message. All three processes exited successfully and wrote
complete tables; **successful exit does not imply DMFT convergence**.

At U=1.5, the first three positive-frequency values of `-Im G` are
`1.81774, 1.53451, 1.33044`; at U=2.5 they are
`1.42815, 0.927608, 0.729623`. Both rise toward low frequency and retain a
central real-frequency peak, consistent with metallic-like solutions rather
than an apparent transition at U=1.5. The U=3.5 run strongly suppresses central
weight but is not converged and cannot be used to assign a phase.

Positive-frequency Matsubara imaginary parts and all sampled retarded imaginary
parts have the causal sign. Matsubara real parts are small (maximum absolute
values approximately `2.73e-7`, `1.10e-6`, and `4.60e-9`, respectively).
The real-axis spectra are not exactly particle-hole symmetric: maximum sampled
`|Im G(w)-Im G(-w)|` is about `0.0872`, `0.0668`, and `0.00249`.
Finite bath size, finite broadening, solver accuracy, and the solver's stopping
criterion limit quantitative interpretation; no precise transition is inferred.

## Artifacts and checks

Run directories, relative to `build_mpi/runs/garcia2004_star_cpu/`:

- `nb7_u1p5_cpu_20261120T181612Z` (input/launcher commit `a5543e6d`)
- `nb7_u2p5_cpu_20261120T181645Z` (input/launcher commit `a5543e6d`)
- `nb7_u3p5_cpu_20261120T181821Z` (input/launcher commit `94e6a266`)

Each contains the input, executable/input checksums and environment in
`manifest.txt`, full launcher and solver logs, and Matsubara/lattice/real-axis
Green-function tables. Generated outputs remain untracked under `build_mpi/`.
A summary can be regenerated with:

```sh
python3 cincuenta/examples/garcia2004_star/summarize_runs.py \
  build_mpi/runs/garcia2004_star_cpu
```

The three small DMRG/MPI smoke checks and five two-iteration MPI preflight
checks passed. The latter compare one- and two-rank equilibrium outputs.
The two analytic DOS/lattice tests passed in the serial CPU build.

The reported postdoc result has not been reproduced here. Next comparisons
should use the actual postdoc input and executable provenance, and distinguish
standard matvec from the GPU BatchedGemm path. A longer U=3.5 run is needed
before interpreting that high-interaction control as converged.
