```@meta
CurrentModule = GenLAProblems
```

# Usage and Quickstart

## Quickstart

```julia
using GenLAProblems

# Build a small system and reduce to REF/RREF.
A = [1 2 3; 2 4 6; 1 1 1]
matrices, pivots, desc = reduce_to_ref(A; gj=true)

# Render a final GE table (requires Python stack).
svg = show_ge_final(matrices, desc, pivots; n_rhs=0)
```

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

## Rendering examples

```julia
using GenLAProblems

# Example system
A = [1 2; 3 4]
b = [5, 6]

# Render system as SVG
pb = ShowGE(A, b; output_dir="/tmp/la")
show_system(pb; var_name="x", fig_scale=1.2)

# Access RHS block (or a single RHS column) from the final GE stack:
rhs = rhs_block(pb)
rhs1 = rhs_block(pb; b_col=1)

# Render backsubstitution cascade
U = [1 2; 0 3]
show_backsubstitution(U, b)

# Render final GE table directly
matrices, pivots, desc = reduce_to_ref(A; gj=true)
show_ge_final(matrices, desc, pivots)
```

Rendering defaults to subdirectories under `/tmp/la` unless you pass
`output_dir` or the compatibility alias `tmp_dir`.

### Inconsistent systems

When a system is inconsistent, `ShowGE` records per-RHS status in
`pb.rhs_status` (e.g., `[:inconsistent, :consistent, ...]`) and sets
`pb.status` to `:inconsistent` or `:mixed`. In the rendered layout, inconsistent
RHS columns are marked with a red **×** in the variable-summary row. The
backsubstitution view is reduced to `0 = rhs` with **No Solution**, and
`show_solution!` returns an empty vector for those RHS columns. The
`solutions(pb)` helper returns only consistent particular solutions.

## Pure Julia vs Python-backed APIs

Pure Julia (no Python needed):
- Matrix generation helpers (`gen_*`, `unit_lower`, `lower`, `ref_matrix`, `rref_matrix`)
- REF/RREF helpers (`reduce_to_ref`, `normal_eq_reduce_to_ref`, `homogeneous_solutions`)
- QR helpers (`gram_schmidt_w`, `normalize_columns`, `gram_schmidt_stable`, `qr_layout`)
- Polynomial helpers (`charpoly`)

Python-backed (requires `LAFigureSpecs`, `matrixlayout`, `jupyter_tikz`, `sympy` via PythonCall):
- Rendering helpers (`show_layout!`, `show_system`, `show_backsubstitution`, `show_solution`, `show_ge_final`)
- SymPy helpers (`sym_*`, `sym_subs_numeric`, `sym_to_julia_*`)

## Adjoint inputs

Most top-level APIs accept `AbstractMatrix` and will work with adjoint wrappers like `A'`.
In-place helpers that mutate matrices (e.g., row swap/eliminate utilities) require mutable
matrices and will not work with `Adjoint`/`Transpose` inputs.
