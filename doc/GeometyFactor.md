# Per-direction `GeometryFactor`

`GeometryFactor` multiplies the geometry contribution to a Hamiltonian link.
The value is evaluated as a DMRG++ input expression at the current simulation
time.  It can therefore be a number, a complex expression, or a time-dependent
expression such as `exp:*:1.0i:!readTable(AversusTime,%t)`.

A factor now belongs to one **geometry term** and one **geometry direction**.
This makes it possible to apply different multipliers to the different bond
directions of the same geometry, and to the same direction in different
geometry terms.

## Input syntax

Use this spelling for a factor associated with geometry term `X` and direction
`Y`:

```text
gtX:dirY:GeometryFactor=<expression>;
```

`X` is the zero-based geometry-term index from `NumberOfTerms`; `Y` is the
zero-based direction index used by that term's `dirY:Connectors` entry.  The
factor applies only to links selected by that geometry term and direction.

For example, a two-term, two-leg ladder can assign a distinct factor to every
term/direction pair:

```text
NumberOfTerms=2;

gt0:GeometryKind="ladder";
gt0:GeometryOptions="ConstantValues";
gt0:LadderLeg=2;
gt0:dir0:Connectors=[1.0];
gt0:dir1:Connectors=[1.0];
gt0:dir0:GeometryFactor=1.10;
gt0:dir1:GeometryFactor=0.90;

gt1:GeometryKind="ladder";
gt1:GeometryOptions="ConstantValues";
gt1:LadderLeg=2;
gt1:dir0:Connectors=[1.0];
gt1:dir1:Connectors=[1.0];
gt1:dir0:GeometryFactor=1.20;
gt1:dir1:GeometryFactor=0.80;
```

The direction numbers are defined by `GeometryKind`; use the same direction
number as the corresponding `Connectors` input.  A chain normally has only
`dir0`, while a ladder has the direction entries required by its geometry.
Do not infer a direction number from a site pair.

## Default and compatibility behavior

An omitted `gtX:dirY:GeometryFactor` has factor `1.0`, so it does not change
that geometry contribution.

The old unqualified form remains supported:

```text
GeometryFactor=<expression>;
```

It is an alias only for `gt0:dir0:GeometryFactor`.  It does **not** apply to
other directions of term 0 or to any direction of another term.  For new input
files, prefer the fully qualified spelling, including when `NumberOfTerms=1`:

```text
gt0:dir0:GeometryFactor=exp:*:1.0i:!readTable(AversusTime,%t);
```

If both forms are present, the explicit
`gt0:dir0:GeometryFactor` value takes precedence.

## Scope

Per-direction factors are read for geometry kinds that use explicit
`dirY:Connectors` directions.  `LongRange` and `General` geometries use their
own geometry data path and do not read `GeometryFactor`.

During finite-system sweeps, DMRG++ normalizes the temporary fringe-site pair
before selecting the direction.  Consequently, a qualified factor is selected
by the physical geometry direction even at the system/environment fringe.
