```@meta
CurrentModule = GenLAProblems
```

# Exact Arithmetic

`GenLAProblems` is designed to generate exact classroom examples rather than
floating-point benchmark matrices.

## What to expect

- Most generators start from small integer data.
- Elimination-oriented examples are constructed so they can be represented with
  exact integers or rationals.
- QR/eigen/SVD examples may still introduce structured rational entries in the
  orthogonal factors, depending on the selected family.

## Typical arithmetic behavior

| Generator family | Typical output style |
| --- | --- |
| `gen_gj_pb`, `gen_inconsistent_gj_pb`, `gen_lu_pb`, `gen_plu_pb`, `gen_ldlt_pb` | Small integers and exact rationals |
| `gen_qr_problem(...; family=:pythagorean)` | Small exact matrices built from structured orthogonal seeds |
| `gen_qr_problem(...; family=:cayley)` | Exact rational orthogonal factors with denser mixing |
| `gen_symmetric_eigenproblem`, `gen_svd_problem` | Exact matrices; orthogonal factors may contain rational structure |

## Displaying exact outputs

Use `LAlatex` when you want notebook-friendly rendering of the generated
matrices and factorizations:

```julia
using GenLAProblems
using LAlatex

A, X, B = gen_gj_pb(3, 4, 3; maxint=3, num_rhs=2)
l_show("A X = B : ", A, X, " = ", B)
```

## SymPy substitution helper

`GenLAProblems` also provides a small helper for substituting numeric values
into symbolic matrices while preserving exact values when possible.

```@docs
GenLAProblems.SymPyHelpers.sym_subs_numeric
```
