# Wave-Function Transformation Acceleration

Wave-function transformation (WFT) supplies an initial vector after the DMRG
basis changes.  The `SolverOptions` entries below select the WFT algorithm; they
do not themselves select a thread count.

## Acceleration selection

`WftOptions` chooses one acceleration mode for each WFT operation:

| Configuration | Selected mode |
| --- | --- |
| one-site DMRG, no WFT acceleration option | patches |
| `SolverOptions` includes `twositedmrg`, with no WFT acceleration option | blocks |
| `SolverOptions` includes `wftAccelPatches` | patches, including with `twositedmrg` |
| `SolverOptions` includes `wftNoAccel` | no acceleration; WFT remains enabled |

The option is spelled `wftAccelPatches`, not `wftAccelPathes`.  There is no
`wftAccelBlocks` SolverOptions token: blocks is the two-site default.  If
several of the above tokens are present, the present implementation applies
assignments in this order: patches, SVD, then no acceleration.  Consequently,
`wftNoAccel` wins over the other two.

`wftAccelSvd` is recognized by `WftOptions`, but it is not a supported usable
mode.  It is not documented in the input checker, its selected operator does
not transform the vector, and its dormant internal implementation reports that
it is not ready.  Do not use it.

`wftNoAccel` is distinct from `nowft`: `wftNoAccel` uses the ordinary WFT
algorithm, whereas `nowft` prevents the normal WFT factory from performing the
transformation.

## `ParallelWftMany`

`ParallelWftMany` is an outer parallel loop over target vectors.  A task
constructs a transformed copy of one target vector and writes that result back;
it does not divide a single target vector among its workers.  It has useful
parallel work only where the caller supplies more than one nonempty target
vector.

The named code section is configured by either of these input labels:

```text
Threads=N;
Threads:ParallelWftMany=N;
```

`Threads=N;` is the default for every named code section.
`Threads:ParallelWftMany=N;` overrides only this outer loop.

The current policy is deliberately as follows:

- With `wftNoAccel` (or the inherited `nowft` condition), `ParallelWftMany`
  uses its resolved named-section thread count.
- With patches or blocks acceleration active, it is forced to one outer worker,
  even when `Threads` or `Threads:ParallelWftMany` requests more workers.

For example, this enables target-vector-level parallelism while making the
inner local WFT serial:

```text
SolverOptions=...,wftNoAccel
Threads=1;
Threads:ParallelWftMany=4;
Threads:WaveFunctionTransfLocal=1;
```

## Inner shared-memory parallelism

The accelerated implementations have their own pthread-backed named sections
when DMRG++ is built with `USE_PTHREADS`:

| WFT mode | Inner named sections |
| --- | --- |
| blocks | `WftAccelBlocks` |
| patches | `BlockDiagWfCtor`, `BlockDiagWfTransform` |
| ordinary, no acceleration | `WaveFunctionTransfLocal` |
| SVD (unfinished) | `WftAccelSvdOne`, `WftAccelSvdTwo` |

For example, tune the default two-site blocks algorithm with:

```text
Threads:WftAccelBlocks=4;
```

or patches with:

```text
Threads:BlockDiagWfCtor=4;
Threads:BlockDiagWfTransform=4;
```

A non-pthread build rejects requests for more than one thread.  The outer
`ParallelWftMany` policy prevents concurrent accelerated WFT calls, but the
inner sections remain independently parallel.  Conversely, a no-acceleration
run can enable both the outer `ParallelWftMany` and the inner
`WaveFunctionTransfLocal` sections.  The code does not coordinate those two
counts, so setting both above one can oversubscribe the available cores.

## Rejected accelerated configurations

DMRG++ rejects a `wftSome` call when its resolved
`Threads:ParallelWftMany` count is greater than one and the selected
accelerated implementation also has an inner count greater than one.  This
uses resolved named-section values, so a global `Threads=N;` request counts
unless a relevant section has its own override.

The conflicting inner sections are:

| WFT mode | Rejected when this section exceeds one |
| --- | --- |
| blocks | `WftAccelBlocks` |
| patches | `BlockDiagWfCtor` or `BlockDiagWfTransform` |

Set the outer section to one when using accelerated inner parallelism:

```text
Threads:ParallelWftMany=1;
Threads:WftAccelBlocks=4;
```

Or keep every applicable accelerated inner section serial when requesting
outer target-vector parallelism:

```text
Threads:ParallelWftMany=4;
Threads:BlockDiagWfCtor=1;
Threads:BlockDiagWfTransform=1;
```

This rejection does not apply to `wftNoAccel`; it does not alter the existing
outer-loop serial fallback for otherwise valid accelerated configurations.
