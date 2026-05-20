```@meta
CurrentModule = GenLAProblems
```

# GenLAProblems

Documentation for [GenLAProblems](https://github.com/ea42gh/GenLAProblems.jl).

## Compatibility matrix

| Component | Minimum | Notes |
| --- | --- | --- |
| Julia | 1.9 | Tested with 1.9–1.10. |
| PythonCall | 0.9 | Required for LAFigureSpecs/matrixlayout interop. |
| LAFigureSpecs | current | Spec generation + algorithms. |
| matrixlayout | current | TeX + SVG rendering. |
| jupyter_tikz | itikz_port | TeX toolchains + SVG. |
| TeX toolchain | TeX Live 2022+ | Needed for SVG rendering. |

## CLI quick checks

Run a smoke render to validate toolchains:

```
python render_smoke.py
```

Set `GENLAPROBLEMS_SMOKE_OUT` to control the output directory.

## Guides

- [Usage and quickstart](usage.md)
- [Python dependencies and rendering](python-deps.md)
- [Exact arithmetic and SymPy helpers](exact-arithmetic.md)

## Generator notes

- `gen_qr_problem` is the canonical QR dispatcher and supports
  `family=:pythagorean`, `:hadamard`, `:cayley`, and `:sparse`.
- `gen_qr_problem_hadamard` is the explicit Hadamard-only QR constructor.
- `gen_svd_problem` supports independent `left_family` / `right_family`
  keywords for the orthogonal factors.
- `W_2_matrix`, `W_3_matrix`, `W_4_matrix`, and `W_matrix` return `(d, A)`
  with diagonal `A' * A`; the corresponding `Q_*` constructors return
  orthogonal matrices.

```@index
```

```@autodocs
Modules = [GenLAProblems]
```
