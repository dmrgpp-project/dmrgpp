# Matched Nb7 ED repeats (2026-11-20)

## Setup

The three `*_ed.ain` inputs repeat the Nb7 CPU DMRG cases at U=1.5, 2.5,
and 3.5, changing only `ImpuritySolver="exactdiag"`, the output root, and a
comment. Bath initialization, target particle sector, frequency grids, fitting
settings, eight-iteration cap, and tolerance remain unchanged. DMRG-specific
keys remain present but do not control ED truncation.

Each ED solve used one MPI rank, one BLAS/OpenMP thread, and the same CPU MPI
executable as the DMRG repeats. The ED path does not distribute frequencies
among MPI ranks. Run input/launcher commit: `a8c27142`.

**No numerical solver correction was made for these repeats.** ED's existing
Matsubara-tail normalization is enabled; DMRG's is disabled. This difference
matters for interpreting the results below.

## Results as produced by the current code

| U | SCF solves | Final reported error | Meets tolerance? | -Im G(i omega_0), rescaled | A(0), unscaled real-axis output | Seconds |
|---|---:|---:|---|---:|---:|---:|
| 1.5 | 3 | 2.20684e-6 | yes | 1.91098 | 2.07358 | 2 |
| 2.5 | 5 | 2.92550e-5 | yes | 1.74884 | 1.55397 | 4 |
| 3.5 | 8 | 3.63563e-3 | **no** | 0.900919 | 0.537638 | 6 |

`omega_0=pi/64`; `A(0)=-Im G(0+i*0.1)/pi`. All three runs exited zero and
wrote complete output tables. The U=3.5 result is not a converged solution.

**These are not clean quantitative ED-versus-DMRG comparisons:** the ED path
has an independently identified factor-four spectral-weight defect, plus a
finite-frequency normalization distortion described next. Absolute ED spectral
amplitudes in the raw comparison plots should not be interpreted physically.

## Factor-four ED spectral-weight defect

The local Green function is for one spin at one impurity site. In
`LanczosPlusPlus/src/Engine/Engine.h`:

1. `spectralFunction` identifies the same-site/same-orbital request as diagonal
   and keeps the two particle/hole types (0 and 1).
2. `getModifiedState` first accumulates the operator at `isite`, then accumulates
   it again at `jsite`. For this request these are the same operator and site,
   with the same sign. The modified vector is therefore doubled. The early
   return for `Tj1Orb.h` does not apply to this Hubbard calculation.
3. `calcSpectral` weights each continued fraction by the squared norm of that
   vector and uses `diagonalFactor=1` for this diagonal request. The resulting
   local Green function retains a factor of four in its spectral weight.

The U=1.5 launcher log directly exposes the doubling on the first solve:

```text
isite=0 type=0 modif=(0.499999642033,0)
jsite=0 type=0 modif=(1.99999856813,0)
isite=0 type=1 modif=(0.500000357967,0)
jsite=0 type=1 modif=(2.00000143187,0)
```

The summed final particle/hole weights are four instead of one. A trapezoidal
integral of the broadened ED `-Im G/pi` over the available range [-4,4] gives
3.93263, 3.92635, and 3.91431 for the three U values. Missing Lorentzian tails
outside the grid account for some weight not included by this finite-window
check; the integrals are consistent with four, not one.

## Additional Matsubara normalization issue

After each Matsubara solve, `ImpuritySolverBase::enforceSpectralSumRule` divides
the Green function by the mean of `-omega*Im G` at the last five frequencies.
The final-iteration diagnostics are:

| U | Reported tail weight | Applied factor | Applied factor times four |
|---|---:|---:|---:|
| 1.5 | 3.68213 | 0.271582 | 1.08633 |
| 2.5 | 3.35301 | 0.298239 | 1.19296 |
| 3.5 | 2.93808 | 0.340358 | 1.36143 |

This cancels the constant four but also rescales the physical Green function
because the last frequency is only about 3.09, not the infinite-frequency
limit. The rescaled Matsubara function feeds the next DMFT self-energy and
bath fit. Thus dividing only the plotted real-axis curve by four would not
repair the underlying self-consistent solution. The final real-axis solve is
not subjected to this normalization and retains its factor-four weight.

These ED findings do not establish the cause of the earlier postdoc plots.
The shared spectral-engine defect and the finite-tail rescaling need separate
regression coverage and correction before using ED as a quantitative reference.

## Artifacts

Relative to `build_mpi/runs/garcia2004_star_ed/`:

- `nb7_u1p5_ed_20261120T205150Z`
- `nb7_u2p5_ed_20261120T205152Z`
- `nb7_u3p5_ed_20261120T205157Z`

Inputs, executable checksums, environment, logs, and output tables are retained
in these directories. Generated plots (no manual amplitude rescaling):

- `build_mpi/plots/nb7-ed-green-functions.svg`
- `build_mpi/plots/nb7-ed-vs-dmrg.svg`

Plots show `-Im G`, **not** `-Im G/pi`, and explicitly warn that ED Matsubara
and real-axis normalization differ. Unconverged solutions have distinct line
styles and labels. The solver-aware summarizer and plotter also continue to
accept the original DMRG runs.
