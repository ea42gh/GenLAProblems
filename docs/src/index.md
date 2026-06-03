```@meta
CurrentModule = GenLAProblems
```

# GenLAProblems

Documentation for [GenLAProblems](https://github.com/ea42gh/GenLAProblems.jl).

## Compatibility matrix

| Component | Minimum | Notes |
| --- | --- | --- |
| Julia | 1.10 | Tested with 1.10–1.12. |
| LAlatex | current | Optional companion package for notebook display. |

## Guides

- [Usage and quickstart](usage.md)
- [Python dependencies and rendering](python-deps.md)
- [Exact arithmetic and SymPy helpers](exact-arithmetic.md)

## Generator notes

- `gen_qr_problem` is the canonical QR dispatcher and supports
  `family=:pythagorean`, `:hadamard`, `:cayley`, and `:sparse`.
- `gen_svd_problem` supports independent `left_family` / `right_family`
  keywords for the orthogonal factors.
- For exact rational SVD factors, `family=:hadamard` additionally requires the
  size to have an integer square root; a Hadamard order alone is not enough.
- `W_2_matrix`, `W_3_matrix`, `W_4_matrix`, and `W_matrix` return `(d, A)`
  with diagonal `A' * A`; the corresponding `Q_*` constructors return
  orthogonal matrices.

For GE reduction workflows, `ShowGE`, and rendered teaching displays, use
`LATeachingSuite` rather than `GenLAProblems`.

```@index
```

```@autodocs
Modules = [GenLAProblems]
```
