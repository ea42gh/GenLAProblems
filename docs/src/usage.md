```@meta
CurrentModule = GenLAProblems
```

# Usage and Quickstart

## Quickstart

```julia
using GenLAProblems

# Build a small exact exercise.
A_ge, X_ge, B_ge = gen_gj_pb(3, 4, 3; maxint=3, num_rhs=2)

# Build an inconsistent GE/GJ system.
A_bad, B_bad = gen_inconsistent_gj_pb(4, 6, 3; maxint=2, num_rhs=1)
```

For Gaussian-elimination workflows, reduction helpers, and teaching displays,
the canonical Julia surface is `LATeachingSuite`.

## Problem generators

### GE/GJ systems

- `gen_gj_pb(m, n, r; ...)` returns a consistent system `A, X, B` with `A * X = B`.
- `gen_inconsistent_gj_pb(m, n, r; ...)` returns `A, B` with each RHS column
  chosen outside `col(A)`, so the system is inconsistent.
- `gen_inconsistent_gj_pb` requires `r < m`, since a full row-rank system
  cannot be inconsistent.

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

- `W_2_matrix`, `W_3_matrix`, `W_4_matrix`, and `W_matrix` return `(d, W)`
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

### Metadata exports

- `__version__` is the package version string for the loaded module.
- `__build__` is build or source-revision metadata when available.

### Basic matrix constructors

- `symbol_vector(s, indices)` returns a vector of `Symbol` names such as `:x_1`, `:x_2`, `:x_3`.
- `symbols_matrix(s, row_indices, col_indices)` returns a matrix of `Symbol` names indexed by the supplied row and column labels.
- `unit_lower(m, n; ...)` and `unit_lower(m; ...)` construct unit lower-triangular matrices.
- `lower(m, n; ...)` and `lower(m; ...)` construct lower-triangular matrices with nonzero diagonal.
- `rref_matrix(m, n, r; ...)` returns a reduced row-echelon matrix and its pivot columns.
- `ref_matrix(m, n, r; ...)` returns a row-echelon matrix and its pivot columns.
- `symmetric_matrix(m; ...)` returns a symmetric integer matrix.
- `skew_symmetric_matrix(m; ...)` returns a skew-symmetric integer matrix.
- `e_i(i, n)` returns the `i`th standard basis vector in `R^n`.
- `i_with_onecol(m, c; ...)` constructs an identity-like elimination matrix with one randomized column.
- `gen_permutation_matrix(row_order)` and `gen_permutation_matrix(n)` construct permutation matrices.
- `gen_full_col_rank_matrix(mc, nc; ...)` builds a small exact matrix with full column rank.
- `ca_projection_matrix(A)` returns the orthogonal projection matrix onto `col(A)` when `A' * A` is invertible.

### Linear-system and factorization helpers

- `gen_gj_matrix(m, n, r; ...)` builds the coefficient matrix used by the GE/GJ problem generators.
- `gen_rhs(A, pivot_cols; ...)` returns a consistent solution/RHS pair `X, B` with `B = A * X`.
- `gen_particular_solution(pivot_cols, n; ...)` returns a particular solution with free variables set to zero.
- `gen_gj_pb(m, n, r; ...)` and `gen_gj_pb(m, n; ...)` generate consistent systems `A, X, B`.
- `gen_inconsistent_gj_pb(m, n, r; ...)` generates inconsistent systems `A, B`.
- `invert_unit_lower(L)` returns the inverse of a unit lower-triangular matrix.
- `gen_inv_pb(n; ...)` returns an invertible matrix together with its exact inverse.
- `gen_lu_pb(m, n, r; ...)` returns `pivot_cols, L, U, A` with `A = L * U`.
- `gen_plu_pb(m, n, r; ..., nswaps=...)` returns `pivot_cols, P, L, U, A` with `A = P * L * U`; it creates dependent rows below rank and then permutes them upward to force up to `min(nswaps, m-r, r-1)` row exchanges. When `r == m`, this construction has no zero rows to work with, so it may return `P = I`.
- `gen_ldlt_pb(m; ...)` returns `L, D, A` with `A = L * D * L'`.

### Orthogonal and QR constructors

- `W_2_matrix`, `W_3_matrix`, `W_4_matrix`, and `W_matrix` return `(d, W)` where `(W // d)` is orthogonal; `W_matrix(n; maxint=...)` forwards the range control to its applicable construction.
- `Q_2_matrix`, `Q_3_matrix`, `Q_4_blocks`, `Q_4_matrix`, and `Q_matrix` return exact orthogonal matrices.
- `sparse_Q_matrix(n; ...)` builds a block-structured exact orthogonal matrix.
- `sparse_W_matrix(n; maxint=...)` returns the denominator-factor representation of `sparse_Q_matrix(n)` and forwards `maxint`.
- `gen_qr_problem(n; ...)` is the canonical QR exercise generator.

### Eigenvalue, Jordan, and SVD constructors

- `gen_eigenproblem(e_vals; ...)` returns `S, Λ, S_inv, A` with prescribed eigenvalues.
- `gen_cx_eigenproblem(evals_no_conj; ...)` returns a real matrix with complex-conjugate eigenstructure encoded in `2 × 2` blocks.
- `gen_symmetric_eigenproblem(e_vals; ...)` returns `Q, Λ, A` with `A = Q * Λ * Q'`.
- `gen_non_diagonalizable_eigenproblem(e_dup, e; ...)` returns a small non-diagonalizable matrix.
- `jordan_block(λ, k)` constructs a single Jordan block.
- `jordan_form(j_blocks)` assembles a block-diagonal Jordan matrix.
- `gen_from_jordan_form(j_blocks; ...)` constructs a matrix similar to a supplied Jordan form.
- `gen_degenerate_matrix(block_descriptions...; ...)` returns `P, J, P_inv, A` with `A = P * J * P_inv`.
- `gen_svd_problem(m, n, σ; ...)` returns `U, Σ, Vt, A` for an exact SVD exercise.

## Adjoint inputs

Most generator APIs accept `AbstractMatrix` and will work with adjoint wrappers
like `A'` when no in-place mutation is required.
