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

```@index
```

```@autodocs
Modules = [GenLAProblems]
```
