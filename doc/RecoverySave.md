# RecoverySave

`RecoverySave` writes restart checkpoints during finite DMRG loops.  It is
intended to make an interrupted calculation restartable by running the same
input again in the same working directory.

## Minimal configuration

Enable checkpoint writing with `RecoverySave=` and enable automatic reading
with `recoveryenableread` in `SolverOptions`:

```text
SolverOptions=twositedmrg,recoveryenableread;
OutputFile="data210";

FiniteLoops=[[@auto, 100, 0],
             [@auto, 200, 0],
             [@auto, 200, 0]
];

string RecoverySave="%l%%1,@keep,@M=200";
```

`RecoverySave` is disabled when the input omits it or sets it to `no`.
Whitespace is not allowed inside its comma-separated value.

The predicate portion controls when a finite-loop checkpoint is written.
It may refer to `%l`, the finite-loop index, `%t`, the targeting time, and
`%ll`, the finite-loop length.  In the example, `%l%%1` writes after every
finite loop.

## RecoverySave options

The options follow the predicate and use `@`:

| Option | Meaning |
| --- | --- |
| `@keep` | Retain recovery files after normal completion. Without it, they are removed when the run finishes normally. |
| `@M=n` | Keep at most `n` recovery slots before overwriting the oldest slot. The default is 10. |

Recovery files are named `Recovery<N><output-file>`.  For example, an
`OutputFile="data210"` calculation writes files such as
`Recovery0data210.hd5` in its working directory.

## Automatic restart

When `SolverOptions` includes `recoveryenableread`, DMRG++ searches the
current working directory for matching recovery files before starting the
solver.  It selects a structurally complete checkpoint and resumes its saved
finite-loop index and sweep position.

Current-format recovery files contain the following metadata under
`/Def/Recovery/`:

- `loopIndex`
- `stepCurrent`
- `generation`
- `complete`

`generation` orders checkpoint writes independently of the ring-slot number.
DMRG++ selects the complete candidate with the greatest generation, so a
wrapped ring uses the last checkpoint written rather than the numerically
largest `Recovery<N>` name.  A file with `generation` but no true `complete`
marker is ignored, which preserves the preceding complete checkpoint if an
interruption occurs during a write.

For compatibility, if no current-format recovery file is available, DMRG++
can select a legacy recovery file with readable loop and step metadata by its
modification time.  It never ranks legacy and current-format checkpoints in
the same ordering.

Automatic recovery uses restart semantics while resolving Ainur `@auto`
finite-loop lengths.  Consequently, unchanged inputs that use both
`recoveryenableread` and `@auto` can resume from a recovery boundary.

## Completion behavior

A checkpoint whose `loopIndex` is at or beyond the resolved finite-loop
schedule is terminal.  For automatic recovery, DMRG++ prints:

```text
RecoverySave: checkpoint already completed; exiting.
```

and exits successfully without constructing the solver or modifying the final
output or retained recovery files.  This makes a repeated batch attempt after
a completed calculation idempotent.

`@keep` preserves files for inspection or explicit use; it does not cause an
unchanged automatic invocation to rerun a completed calculation.

## Explicit restart is separate

`SolverOptions=restart` with `RestartFilename=` is an explicit restart
request.  It takes precedence over automatic recovery-file discovery, even if
matching `Recovery<N>` files are present in the directory.  Automatic
terminal-checkpoint no-op behavior applies only to `recoveryenableread`'s
automatic path.

## Operational guidance

- Use a distinct working directory for each logical calculation. Matching
  recovery files from an earlier calculation can otherwise be considered.
- Set `@M` large enough that the checkpoint interval and ring size suit the
  expected interruption and walltime policy.
- Treat recovery files as a set. Do not manually rename or partially copy a
  checkpoint if its associated serialized stacks are required.
- `SolverOptions=noClobber` is not a RecoverySave completion check. It uses
  the redirected text log's normal-completion marker; it does not inspect HDF5
  recovery metadata.

The regression fixture is `TestSuite/inputs/input210.ain`. Its tests are
registered as `recovery_auto_resume`, `recovery_terminal_rerun`, and
`recovery_ring_wrap`.
