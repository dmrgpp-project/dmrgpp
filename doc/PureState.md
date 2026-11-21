# Pure-state initial states

DMRG++ can initialize the infinite-system growth from a specified product
state, rather than obtaining the initial vector by diagonalizing the
superblock Hamiltonian.  The input is a `gs` expression of the form
`PureState[...]`.

This feature is useful when a calculation must start from a known occupation
configuration and subsequently use the normal finite-sweep or restart
workflow.

## Input syntax

Add one `gs` line to the Ainur input:

```text
gs=PureState[i0,i1,...,iN-1];
```

There must be exactly one index for every physical site, in physical-site
order.  Each index selects the local basis state for its site; it must be a
non-negative decimal integer valid for that site's local Hilbert space.

For example, `TestSuite/inputs/input160.ain` specifies a six-site product
state:

```text
TotalNumberOfSites=6;
SolverOptions="twositedmrg,geometryallinsystem";
# ... model and DMRG parameters ...
gs=PureState[0,0,2,3,0,1];
```

The grammar is deliberately strict:

- the expression must begin with `PureState[` and end with `]`;
- the comma-separated list must contain at least one index;
- every index must contain decimal digits only; do not include signs or
  whitespace inside the list; and
- the number of indices must equal `TotalNumberOfSites`.

An invalid expression is rejected during input processing.  An index outside
the local Hilbert-space range is a user error and must not be supplied.

## Symmetries and growth geometry

DMRG++ derives the quantum-number adjustments required by the product state
from the selected local indices.  Therefore, do **not** set
`AdjustQuantumNumbers` together with `gs=PureState[...]`; the input parser
rejects that combination.

Both supported infinite-growth layouts use the same `gs` syntax:

- With `geometryallinsystem` in `SolverOptions`, growth is all in the system.
  The product state must contain at least three sites.
- Without `geometryallinsystem`, growth is symmetric between system and
  environment.  The product state must contain an even number of sites and
  at least four sites.

The implementation uses the model's local quantum numbers to construct the
successive symmetry sectors for the chosen layout.  The site indices are not
occupation numbers in a model-independent sense: consult the model's local
basis convention when selecting them.

## Calculation flow

During the infinite-system stage, the solver transports the requested product
state through each basis truncation and uses the resulting product vector as
the target state.  It does not diagonalize the superblock Hamiltonian in that
stage.  The normal diagonalization path is used again during finite sweeps.

A run initialized this way writes its checkpoint in the usual manner.  It can
be restarted with the normal `restart` and `RestartFilename` settings; see
`TestSuite/inputs/input161.ain` for an example that restarts the
`input160.ain` output before time evolution.

## Symmetric-growth example

The following is the corresponding symmetric-growth form from
`TestSuite/inputs/input162.ain`:

```text
TotalNumberOfSites=6;
SolverOptions="twositedmrg";
# ... model and DMRG parameters ...
gs=PureState[3,2,0,0,0,1];
```

## Diagnostics

The parser reports malformed expressions as `Expected PureState[...] for gs`
and invalid list members as `PureState indices must be non-negative integers`.
During growth, `PureState` progress messages identify the selected local index
for the system and environment sites and report the norms of the transported
vectors.  These messages can help locate an input whose selected local state
is incompatible with the model or symmetry sector.
