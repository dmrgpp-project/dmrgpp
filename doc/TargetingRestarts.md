# Targeting restarts

## Scope

This note defines the restart-format and compatibility contract implemented by
DMRG++.  Restart files identify their source targeting implementation, and the
reader rejects every transition outside the explicit compatibility relation.
Restarting older files without that identity is unsupported.

Terminology used below:

- **source targeting** is the implementation that wrote the restart file;
- **destination targeting** is the implementation selected for the new run;
- **ground-state targeting** means the implicit default when no explicit
  targeting option is present, as well as explicit `GroundStateTargeting`.

## How a restart is selected and performed

`TargetSelector` resolves the destination targeting implementation from the
new run's `SolverOptions`.  Before that implementation is constructed or any
restart state is read, `DmrgSolver::main` opens the restart file and validates
the persisted source ID against the resolved destination ID
(`dmrg/Engine/DmrgSolver.h:234-244`, `TargetSelector.h:93-170`).  No explicit
targeting option resolves to the canonical `GroundState` ID.

Only after that check does `psi.read(io, "FinalPsi")` dispatch to the
destination reader.  An incompatible pair therefore fails before checkpoint
or target-vector state is consumed.

Every final and recovery writer creates `FinalPsi/Provenance/TargetingId`.
It is a mandatory canonical implementation ID, independent of the serialized
target-vector representation (`DmrgSolver.h:478-482`, `Recovery.h:258-264`).
Every targeting writer also stores common wavefunction state using
`TargetingCommon::write`: `TargetCentralSite` plus `PSI`
(`TargetingCommon.h:227-242`).  Targetings needing target vectors, stages, or
time additionally store `FinalPsi/TimeSerializer`.  The serializer records
its writer-supplied name, current time and time step, stages, and target
vectors (`TimeSerializer.h:99-117,203-227`).

The generic time-aware reader, `TargetingCommon::readGSandNGSTs`, always reads
the common state and then behaves as follows (`TargetingCommon.h:284-374`):

1. A missing `TimeSerializer` is accepted silently only for an approved
   `GroundState` source.  Its additional destination state remains at newly
   constructed defaults.  For every other source, absence is an error.
2. It maps source target vectors through `RestartMappingTvs`; invalid
   non-negative mappings fail, while negative mappings leave the destination
   vector untouched.
3. It restores stages only when `RestartMapStages` is enabled, requiring the
   source and destination stage counts to agree.
4. It restores global time and time-step only when the stored serializer name
   equals the destination reader's name (or, for old serializers with no name,
   when the target-vector counts agree).  Otherwise it resets both to zero
   (`TargetingCommon.h:368-373,785-794`).

The source/destination policy is separate from this vector-format mapping and
is enforced before this reader executes.

## Compatibility relation

`TargetingBase::restartIsCompatible` (`TargetingBase.h:356-369`) accepts:

1. a self-restart with the same canonical ID;
2. `GroundState -> <any destination>`;
3. `CorrectionVector -> RixsStatic`; and
4. `RixsStatic -> RixsDynamic`.

All other source/destination pairs fail with both canonical IDs in the error
message.  A ground-state restart contains the common `PSI` state but no
targeting-specific serialized state, so it is admissible input for every
destination.  This compatibility decision remains separate from checks for
the WFT, operators, model, and targeting-specific input required by that
destination.

The RIXS adapters justify the two specialised edges: static RIXS requires four
serialized target vectors and restores vectors 1--3, while dynamic RIXS
requires six (`TargetingRixsStatic.h:214-226`,
`TargetingRixsDynamic.h:341-353`).

## Canonical IDs and legacy files

The persisted values are implementation IDs, not user option spellings:
`GroundState`, `TimeStep`, `Dynamic`, `CorrectionVector`, `Correction`,
`Metts`, `RixsStatic`, `RixsDynamic`, `Chebyshev`, `Expression`, and
`CVEvolution`.  `TimeStepTargeting` and `TargetingAncilla` resolve to the same
`TimeStep` ID (`TargetSelector.h:145-170`).

There is no compatibility fallback for a file that lacks
`FinalPsi/Provenance/TargetingId`.  Such a file fails before destination state
is constructed.  The reader does not infer an identity from `SolverOptions`,
`InputBase64Encoded`, or `TimeSerializer/Name`.

`TargetingCVEvolution::write` now writes the `CVEvolution` serializer name,
matching its read path (`TargetingCVEvolution.h:226-242`).  Its self-restarts
therefore no longer reset global time solely because of the previous
`TimeStep` name mismatch.

## Regression coverage

The `input9100 -> input9102` producer/consumer chain exercises one permitted
`GroundState -> Expression` transition.  `input9100_targeting_id` inspects the
written HDF5 field and requires `GroundState`.  A build-tree input derived from
`input9103.inp` attempts `Expression -> GroundState`; the
`incompatible_targeting_restart` test requires the rejection diagnostic before
the destination reader runs (`dmrg/CMakeLists.txt:411-450`).
