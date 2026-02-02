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
svg = show_ge_final(matrices, desc, pivots; Nrhs=0)
```

## Pure Julia vs Python-backed APIs

Pure Julia (no Python needed):
- Matrix generation helpers (`gen_*`, `unit_lower`, `lower`, `ref_matrix`, `rref_matrix`)
- REF/RREF helpers (`reduce_to_ref`, `normal_eq_reduce_to_ref`, `homogeneous_solutions`)
- QR helpers (`gram_schmidt_w`, `normalize_columns`, `gram_schmidt_stable`, `qr_layout`)
- Polynomial helpers (`charpoly`)

Python-backed (requires `la_figures`, `matrixlayout`, `jupyter_tikz`, `sympy` via PythonCall):
- Rendering helpers (`show_layout!`, `show_system`, `show_backsubstitution`, `show_solution`, `show_ge_final`)
- SymPy helpers (`sym_*`, `sym_subs_numeric`, `sym_to_julia_*`)

## Adjoint inputs

Most top-level APIs accept `AbstractMatrix` and will work with adjoint wrappers like `A'`.
In-place helpers that mutate matrices (e.g., row swap/eliminate utilities) require mutable
matrices and will not work with `Adjoint`/`Transpose` inputs.
