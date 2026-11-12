# TimeEvolve

`TimeEvolve` is a global operator for `TargetingExpression`.  It evolves a
P-vector in real time while DMRG++ sweeps.  Enable expression targeting and
complex arithmetic:

```ainur
SolverOptions="twositedmrg,TargetingExpression,useComplex";
```

A declaration has the following form:

```ainur
string P1="TimeEvolve{tau=0.1,steps=5,advanceEach=6}*|P0>";
```

It means that `P1` is the time-evolved target constructed from `P0`.  The
source may also be `|gs>`.  The expression must contain exactly one `*`:
`TimeEvolve{...}*|P...>`.

## Example: `input124.ain`

`TestSuite/inputs/input124.ain` restarts a calculation, reuses the saved
P-vector as its source, and evolves it:

```ainur
SolverOptions="restart,twositedmrg,TargetingExpression,useComplex";
RestartFilename="data123";

string P0=|P1>;
string P1="TimeEvolve{tau=0.1,steps=5,advanceEach=6}*|P0>";
GsWeight=0.2;
```

Here `P0=|P1>` copies P-vector 1 from the restart into the current first
P-vector slot.  `TimeEvolve` then uses that copied state as its source and
places the evolving result in `P1`.

## Parameters

Parameters are comma-separated `key=value` pairs inside braces.  The current
keys are:

| Key | Default | Meaning |
| --- | --- | --- |
| `tau` | `0.1` | Time increment applied when the evolution advances. |
| `steps` | `3` | Number of equally spaced target vectors used for one time interval. It must be greater than one. Their offsets span 0 through `tau`. |
| `advanceEach` | number of sites minus 2 | Number of eligible evolution calls between time advances. A positive advance occurs only at a sweep border. Set `0` to disable advancement. `AdvanceEach` and `advanceeach` are accepted aliases. |
| `startLoop` | `0` | Zero-based finite-loop index before which this evolution does not advance. Its advancement counter begins at the first eligible call in `startLoop`; source availability and `depends` can delay it further. |
| `algorithm` | `Krylov` | Propagator: `Krylov`, `RungeKutta`, `SuzukiTrotter`, or `Chebyshev:a:b`. The Chebyshev form supplies its two transform parameters. |
| `EnergyForExp` | current ground-state energy | Energy origin used for the exponential. `energyforexp` is also accepted. |
| `disposition` | `0` | With `disposition=2`, the source P-vector is destroyed after the first successful time advance. Other values retain it. |
| `depends` | none | Colon-separated prerequisite P-vectors, for example `depends=P0:P3`. The evolution does not advance until every prerequisite is complete. |

For example, the declaration in `input124.ain` makes five target vectors at
times `t`, `t + 0.025`, `t + 0.05`, `t + 0.075`, and `t + 0.1`; it advances
the base time by `0.1` after every six eligible calls when a sweep reaches a
border.

## Scheduling and equal-time measurements

`startLoop` delays an evolution until the requested finite-loop index. Finite
loops are numbered from zero, so `startLoop=2` allows the first clock tick
only in the third finite loop or later:

```ainur
string P1="TimeEvolve{tau=0.1,steps=5,advanceEach=6,startLoop=2}*|P0>";
```

If the requested loop is beyond the configured finite-loop schedule, the base
time never advances during that run.

Every `TimeEvolve` declaration owns its own clock. An evolution begins only
when its source vector is available, so expressions whose sources become
available at different positions in a sweep can otherwise advance at different
times.  DMRG++ rejects a braket involving two time-evolved P-vectors if their
times differ.

When such P-vectors must be measured at equal physical time, give both
evolutions the same explicit prerequisite.  For example, if `P0` is produced
locally, use:

```ainur
string P1="TimeEvolve{tau=0.1,steps=5,advanceEach=6,depends=P0}*|gs>";
string P2="TimeEvolve{tau=0.1,steps=5,advanceEach=6,depends=P0}*|P0>";
```

This creates a shared start barrier: neither evolution advances before `P0` is
complete.  Do not compensate for a delayed source by choosing a different
`advanceEach`; that depends on sweep layout and is not a general synchronization
mechanism.

## Measurements

The primary P-vector name denotes the state at the base time.  For a
`TimeEvolve` target, `|P1.last>` denotes the final auxiliary Krylov slot for
that time interval.  It can be used in a braket, for example:

```text
<P1.last|P1>
```

The regression fixture `TestSuite/inputs/input9103.inp` uses this form along
with an in-situ correlator.

## Limits

`TimeEvolve` is only recognized as a global operator in a
`TargetingExpression`.  It cannot carry a site suffix.  Unrecognized parameter
keys are errors.
