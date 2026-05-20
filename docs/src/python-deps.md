```@meta
CurrentModule = GenLAProblems
```

# Python Dependencies and Rendering

`GenLAProblems` is generator-focused and does not require the Python render
stack for its canonical public API.

If you want GE workflows, rendered teaching figures, or Python-backed display
helpers, use `LATeachingSuite`. That package owns the Julia-facing bridge to:

- `LAFigureSpecs`
- `matrixlayout`
- `jupyter_tikz`
- `sympy`

In other words:

- `GenLAProblems` generates exact example data.
- `LATeachingSuite` handles the cross-language workflow/rendering layer.
