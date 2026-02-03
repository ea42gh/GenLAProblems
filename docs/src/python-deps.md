```@meta
CurrentModule = GenLAProblems
```

# Python Dependencies and Rendering

GenLAProblems uses PythonCall to access the rendering stack and SymPy.
The following packages are expected in the active Python environment:

- `la_figures`
- `matrixlayout`
- `jupyter_tikz`
- `sympy`

If these are missing, rendering helpers like `show_system`, `show_layout!`,
`show_backsubstitution`, `show_solution`, and `show_ge_final` will error.
Make sure PythonCall points to a Python with the required packages (for example,
set `JULIA_PYTHONCALL_EXE`).

## Rendering options (`render_opts`)

High-level rendering helpers accept a `render_opts` keyword. This dictionary is passed
through to the underlying renderers (`la_figures` → `matrixlayout` → `jupyter_tikz`)
and can include options such as:

- `crop`
- `padding`
- `toolchain_name`
- `frame`
- `exact_bbox`

## Troubleshooting

**PythonCall can't find required packages**
- Ensure PythonCall points to the Python where packages are installed:
  `ENV["JULIA_PYTHONCALL_EXE"] = "/path/to/python"`
- Install required packages in that Python:
  `python -m pip install sympy`
  (and install your local `la_figures`, `matrixlayout`, `jupyter_tikz` packages)

**`ModuleNotFoundError: No module named 'la_figures'`**
- Verify `PYTHONPATH` includes the local repos or install the packages into the Python env.
- In tests, `PYTHONPATH` is set to include `0_ITIKZ/la_figures` and `0_ITIKZ/matrixlayout`.

**`ModuleNotFoundError: No module named 'matrixlayout'`**
- Same fix as above; ensure the module is on `PYTHONPATH` or installed.

**`Could not import Python module sympy`**
- Install sympy into the PythonCall interpreter or set `JULIA_PYTHONCALL_EXE` to a Python
  environment that already has sympy installed.
