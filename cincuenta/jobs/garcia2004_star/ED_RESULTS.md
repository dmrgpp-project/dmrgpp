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

## Why the existing regressions missed this

An audit and rerun of the current tests found:

- The two `u0_gImpEqGcluster_*_check` tests compare each solver's impurity
  Matsubara output with its own lattice output, **not ED directly with DMRG**.
  The ED run uses beta 20 and 200 positive Matsubaras; DMRG uses beta 10 and
  six positive Matsubaras. Both perform only one DMFT iteration.
- `DmftSolver` applies ED's sum-rule rescaling before writing the Matsubara
  file tested by the comparator. Thus the test sees a compensated result,
  not the raw continued-fraction weight. The real-axis output is not compared.
- In the local rerun, ED logged `SpectralWeight=3.998990` and
  `NormalizationFactor=0.250063`. Its post-normalization maximum difference
  from the lattice was `0.010110`, passing tolerance `0.015`. Undoing the
  printed normalization factor gives a raw difference of approximately
  `2.77238`. The test's highest Matsubara frequency is about 62.7, so finite-tail
  bias is much smaller than in the Nb7 repeats, whose highest frequency is 3.09.
- The DMRG run completed, but its comparison **failed** locally: maximum
  difference `0.109934` exceeds tolerance `0.05`. Both its run and check are
  labeled `Nightly`; `.github/workflows/ci.yml` runs `ctest -LE Nightly`,
  excluding them from fast CI. This failure was not diagnosed or fixed here.
- Both `spectralWeight` unit tests passed. They exercise synthetic `1/(i omega)`
  and `0.7/(i omega)` vectors, not ED or DMRG solver output or real-axis weights.
- The MPI preflights compare DMRG with itself across rank counts, not with ED.
  The simple `cincuenta0`/`cincuenta1` run tests do not assert spectral amplitudes.

The local rerun log is `build_mpi/normalization-test-audit.log`; the U=0 outputs
are under `build_mpi/cincuenta/src/u0ExactDiag/` and `u0Dmrg/`.

In the current branch's ancestry, `7755170a` (2026-06-19) introduced the
paired U=0 checks and the present tail-normalization/DMRG-opt-out approach.
Related development-branch history further explains the earlier resolution:
`2e821f6c` (2026-05-29) introduced a hardcoded quarter factor; `93b90caa`
(2026-06-03) stopped scaling DMRG. The message for `9a7fcf76` (2026-06-08)
explicitly reports ED's fourfold mismatch and describes sum-rule rescaling as
its fix. `de391243` (committed 2026-06-24) replaced the hardcoded factor with
finite-tail normalization and explicitly opted DMRG out. These latter four
commits are available on local refs but are not ancestors of the current HEAD;
they are development-history evidence, not its literal integration sequence.
The current spectral engine still contains the same-site duplication.

Coverage needed to prevent recurrence: check raw particle/hole weights sum to
one before any correction; compare a known noninteracting/atomic Green function
on both frequency axes; and compare ED and DMRG on an identical fixed bath and
frequency grid, independently of DMFT fitting and output renormalization.
