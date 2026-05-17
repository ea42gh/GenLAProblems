```@meta
CurrentModule = GenLAProblems
```

# Python Dependencies and Rendering

GenLAProblems uses PythonCall to access the rendering stack and SymPy.
The following packages are expected in the active Python environment:

- `LAFigureSpecs`
- `matrixlayout`
- `jupyter_tikz`
- `sympy`

If these are missing, rendering helpers like `show_system`, `show_layout!`,
`show_backsubstitution`, `show_solution`, and `show_ge_final` will error.
Make sure PythonCall points to a Python with the required packages (for example,
set `JULIA_PYTHONCALL_EXE`).

You can explicitly initialize the Python stack with:

```julia
GenLAProblems.ensure_pythoncall!()
```

### Smoke render helper

For a quick toolchain sanity check from Python:

```
python render_smoke.py
```

Set `GENLAPROBLEMS_SMOKE_OUT` to control the output directory.

## Rendering options (`render_opts`)

High-level rendering helpers accept a `render_opts` keyword. This dictionary is passed
through to the underlying renderers (`LAFigureSpecs` → `matrixlayout` → `jupyter_tikz`)
and can include options such as:

- `crop`
- `padding` (tuple order is `(left, right, top, bottom)` in SVG units; use `mm_to_px` to convert mm)
- `toolchain_name`
- `frame`
- `exact_bbox`

`output_dir` is the canonical artifact directory (for TeX/PDF/SVG); `tmp_dir` is a legacy alias.

## Troubleshooting

**PythonCall can't find required packages**
- Ensure PythonCall points to the Python where packages are installed:
  `ENV["JULIA_PYTHONCALL_EXE"] = "/path/to/python"`
- Install required packages in that Python:
  `python -m pip install sympy`
  (and install your local `LAFigureSpecs`, `matrixlayout`, `jupyter_tikz` packages)

**`ModuleNotFoundError: No module named 'LAFigureSpecs'`**
- Verify `PYTHONPATH` includes the local repos or install the packages into the Python env.
- In tests, `PYTHONPATH` is set to include `0_ITIKZ/LAFigureSpecs` and `0_ITIKZ/matrixlayout`.

**`ModuleNotFoundError: No module named 'matrixlayout'`**
- Same fix as above; ensure the module is on `PYTHONPATH` or installed.

**`Could not import Python module sympy`**
- Install sympy into the PythonCall interpreter or set `JULIA_PYTHONCALL_EXE` to a Python
  environment that already has sympy installed.

**`PythonCall disabled` or `unavailable during precompile`**
- Unset `GENLAPROBLEMS_DISABLE_PYTHONCALL` or set `JULIA_PYTHONCALL_EXE` to a working Python.
- Call `GenLAProblems.ensure_pythoncall!()` from a runtime context (not precompile).
