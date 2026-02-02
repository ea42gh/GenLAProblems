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
