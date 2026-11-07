# Finite-loop third entry

In an Ainur `FiniteLoops` row, the third entry controls per-sweep behavior:

```text
FiniteLoops = [
    [@auto, 100, 0],
    [@auto, 200, @save]
];
```

It accepts either a non-negative decimal bit field or one or more named
options separated by `:`.  For example, `@save` is equivalent to `1`, and
`@save:@randomguess` is equivalent to `9`.

## Options

| Name | Bit | Value | Effect |
| --- | ---: | ---: | --- |
| `@save` | 0 | 1 | Save the finite-sweep state to the configured output file. |
| `@onlyfastwft` | 1 | 2 | Use the wave-function-transformed vector and skip diagonalization. The previous energy is retained. |
| `@onlyslowwft` | 2 | 4 | Use the wave-function-transformed vector and skip diagonalization. Its energy is evaluated as the expectation value of the Hamiltonian. |
| `@randomguess` | 3 | 8 | Use a random vector rather than a wave-function-transformed initial guess. |
| `@multisitepush` | 4 | 16 | Send the finite-sweep serializer to the active targeting strategy's `multiSitePush` operation. |
| `@onesitetruncation` | 5 | 32 | Enable the model's one-site-truncation update during this finite loop. |

Named options are lowercase.  A third entry may use named options, a numeric
bit field, and per-loop truncation settings, but it may not combine a numeric
bit field with named options.

## Constraints

- At most one of `@onlyfastwft`, `@onlyslowwft`, and `@randomguess` may be
  set.  This also applies to their numeric bits 1, 2, and 3.
- `@onesitetruncation` requires `@randomguess` (numeric value `40` when both
  bits are represented numerically).
- `@save` can be combined with any option that is otherwise valid.

## Per-loop truncation settings

The third entry can also override truncation-control values for that loop:

| Setting | Meaning |
| --- | --- |
| `tol=<value>` or `tolerance=<value>` | Set the truncation tolerance. |
| `mmin=<value>` | Set the minimum number of kept states. |

Separate each token with `:`.  For example:

```text
FiniteLoops = [
    [@auto, 400, @save:tol=1e-10:mmin=200]
];
```

## Examples

```text
# Ordinary finite sweep; do not save its state.
[@auto, 200, 0]

# Ordinary finite sweep; save its state.
[@auto, 200, @save]

# Save and use a random initial vector before diagonalization.
[@auto, 200, @save:@randomguess]

# Save and use the fast WFT-only path.
[@auto, 200, 3]

# One-site truncation requires a random guess.
[@auto, 200, @randomguess:@onesitetruncation]
```