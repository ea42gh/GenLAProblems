```@meta
CurrentModule = GenLAProblems
```

# GenLAProblems

Documentation for [GenLAProblems](https://github.com/ea42gh/GenLAProblems.jl).

## Rendering options

High-level rendering helpers (e.g., `show_layout!`, `show_system`, `show_backsubstitution!`,
`show_solution!`, and `matrixlayout_ge`) accept a `render_opts` keyword. This dictionary
is passed through to the underlying Python renderers (`la_figures` → `matrixlayout` → `jupyter_tikz`)
and can include options like `crop`, `padding`, `toolchain_name`, `frame`, and `exact_bbox`.

```@index
```

```@autodocs
Modules = [GenLAProblems]
```
