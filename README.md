# GenLAProblems

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ea42gh.github.io/GenLAProblems.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ea42gh.github.io/GenLAProblems.jl/dev/)
[![Build Status](https://github.com/ea42gh/GenLAProblems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ea42gh/GenLAProblems.jl/actions/workflows/CI.yml?query=branch%3Amain)

GenLAProblems depends on `LAlatex` for LaTeX rendering utilities and display helpers.

## Overview

GenLAProblems provides linear algebra problem generators and GE/GJ helpers.

## Display backend

GE/GJ visualizations are rendered via the Python `la_figures` + `matrixlayout` stack
through `PythonCall`. Ensure those Python packages are available in the active
Python environment.

## SymPy helpers (optional)

When mixing Julia arrays with SymPy objects, PythonCall does not auto-convert
types the way PyCall does. For convenience, GenLAProblems provides a small
optional helper module:

```julia
using GenLAProblems.SymPyHelpers

P = circular_shift_matrix(length(v))
Pv = sym_mul(P, v)
@show sym_vec_zero(circular_shift(v) .- Pv)
```

`sympy_subs_numeric` is useful for substitution: it returns a SymPy matrix if
symbols remain, and a Julia numeric array once fully numeric.

Import this submodule only when needed; it is not re-exported by default.

## Local development

If you are working from this monorepo, install the local `LAlatex` checkout in the
`GenLAProblems` environment:

```julia
import Pkg
Pkg.activate("/home/lab/NOTEBOOKS/LSHOW/GenLAProblems")
Pkg.develop(path="/home/lab/NOTEBOOKS/LSHOW/LAlatex")
```
