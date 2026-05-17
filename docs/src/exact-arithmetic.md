```@meta
CurrentModule = GenLAProblems
```

# Exact Arithmetic and SymPy Helpers

## Exact QR helpers

`gram_schmidt_stable` supports exact inputs (e.g., `Rational`, `Complex{Rational}`), but the
exact path returns matrices with symbolic square roots (via Symbolics). This means results are
`Matrix{Any}` and may contain `Symbolics.sqrt(...)` terms. If you need purely numeric results,
use floating inputs or post-process the symbolic output.

`gram_schmidt_w` produces exact, rational-valued orthogonalized columns for integer inputs and
`normalize_columns` introduces symbolic square roots as needed. These helpers are intended for
pedagogical output and may not be suitable for numeric downstream computations without
additional conversion.

## Exact vs numeric (quick comparison)

| Topic | Exact inputs (`Rational`, `Complex{Rational}`) | Numeric inputs (`Float64`, `ComplexF64`) |
| --- | --- | --- |
| `gram_schmidt_stable` | Returns symbolic square roots; `Matrix{Any}` | Returns numeric `Q, R` |
| `gram_schmidt_w` | Rational orthogonal columns | Not intended for floats |
| `normalize_columns` | Uses `Symbolics.sqrt` when needed | Returns numeric Q |
| Downstream math | May need symbolic handling | Direct numeric workflows |

## SymPy helpers

`sym_mat`/`sym_vec` convert Julia arrays to SymPy matrices and preserve exact rationals
(`Rational`) and complex rationals (`Complex{Rational}`) as SymPy rationals and `I`.
If you pass mixed-type `Matrix{Any}` inputs, SymPy conversion can still be sensitive to
PythonCall’s ndarray handling; prefer uniform numeric or symbolic arrays.

`sym_subs_numeric` returns a Julia array only when all symbols are substituted; otherwise it
returns a SymPy matrix. Use `sym_to_julia_mat` or `sym_to_julia_vec` to force conversion when
appropriate.

```@docs
sym_subs_numeric
```

## Known limitations

- Mixed-type `Matrix{Any}` inputs can trigger PythonCall conversion issues in some environments.
- Some rendering helpers require the full Python stack; without it they will throw errors.
- Symbolic outputs from exact QR are not guaranteed to simplify unless you apply Symbolics/SymPy simplification.
