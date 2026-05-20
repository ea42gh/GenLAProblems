```@meta
CurrentModule = GenLAProblems
```

# Usage and Quickstart

## Quickstart

```julia
using GenLAProblems

# Build a small exact exercise.
A_ge, X_ge, B_ge = gen_gj_pb(3, 4, 3; maxint=3, num_rhs=2)
```

For Gaussian-elimination workflows, reduction helpers, and teaching displays,
the canonical Julia surface is `LATeachingSuite`.

## Problem generators

### QR families

`gen_qr_problem` is the canonical QR exercise generator:

```julia
using GenLAProblems

A_pyth = gen_qr_problem(3; family=:pythagorean, maxint=2)
A_had  = gen_qr_problem(4; family=:hadamard, maxint=2)
A_cay  = gen_qr_problem(5; family=:cayley, maxint=2)
A_sp   = gen_qr_problem((2, 2); family=:sparse, maxint=2)
```

Family notes:
- `:pythagorean` uses the small exact `W_k` / `Q_k` families.
- `:hadamard` uses the Hadamard-based QR constructor.
- `:cayley` uses a dense Cayley-transform orthogonal seed.
- `:sparse` uses `sparse_Q_matrix`; when you pass a tuple like `(2, 2)`, it is
  interpreted as block sizes.

### SVD orthogonal-factor families

`gen_svd_problem` now lets you choose the orthogonal-family construction for
the left and right factors independently:

- `m` and `n` may be plain integer dimensions such as `4`, or block-size
  partitions such as `[2, 1]` or `(2, 2)`.
- When you pass partitions, those block sizes are used by the sparse/block
  orthogonal families to build more structured exact factors.

```julia
U, Σ, Vt, A = gen_svd_problem(
    4, 4, [3, 1, 0, 0];
    left_family=:hadamard,
    right_family=:cayley,
    maxint=2,
)
```

The default remains:

```julia
U, Σ, Vt, A = gen_svd_problem([2, 1], [2, 1], [3, 1, 0]; maxint=2)
```

which uses `left_family=:sparse` and `right_family=:sparse`.

For exact rational orthogonal factors, `family=:hadamard` has an extra
constraint: the dimension must both support a Hadamard matrix and have an
integer square root. For example, `4` works, but `12` does not.

### `W_*` and `Q_*` constructors

- `W_2_matrix`, `W_3_matrix`, `W_4_matrix`, and `W_matrix` return `(d, A)`
  where `A' * A` is diagonal.
- `Q_2_matrix`, `Q_3_matrix`, and `Q_4_matrix` return orthogonal matrices.
- `W_matrix` now follows the same two-value contract in both the specialized and
  general cases.

### Eigenvalue input shape

`gen_eigenproblem` and `gen_symmetric_eigenproblem` accept:
- a standard vector like `[3, -1, 2]`
- a row matrix like `[3 -1 2]`
- a column matrix like `[3, -1, 2;;]`

General matrix-shaped eigenvalue input is rejected.

## Display and workflow packages

`GenLAProblems` only owns the generation layer.

- Import `LAlatex` separately when you want to display generated matrices or
  factorizations in notebooks.
- Import `LATeachingSuite` when you want Gaussian-elimination workflows,
  reduction helpers, `ShowGE`, or rendered teaching displays.

## Public surface

Generator-oriented APIs in `GenLAProblems`:
- Matrix generation helpers (`gen_*`, `unit_lower`, `lower`, `ref_matrix`, `rref_matrix`)
- Matrix-construction helpers used by those generators (`W_*`, `Q_*`,
  `sparse_Q_matrix`, `gen_permutation_matrix`, `gen_full_col_rank_matrix`)

## Adjoint inputs

Most generator APIs accept `AbstractMatrix` and will work with adjoint wrappers
like `A'` when no in-place mutation is required.
