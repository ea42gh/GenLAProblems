# GenLAProblems

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ea42gh.github.io/GenLAProblems.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ea42gh.github.io/GenLAProblems.jl/dev/)
[![Build Status](https://github.com/ea42gh/GenLAProblems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ea42gh/GenLAProblems.jl/actions/workflows/CI.yml?query=branch%3Amain)

GenLAProblems depends on `LAlatex` for LaTeX rendering utilities and display helpers.

## Migration status

GenLAProblems is the new home for linear algebra problem generation and GE/GJ helpers.
The legacy `GenLinAlgProblems` package is treated as read-only during the split.

## Display backend transition

Current GE/GJ rendering uses the Python `itikz`/`itikz.nicematrix` stack via `PythonCall`.
We plan to migrate these visualizations to `la_figures` in a later cleanup pass.
Until then, keep `itikz` available in the active Python environment.

## Local development

If you are working from this monorepo, install the local `LAlatex` checkout in the
`GenLAProblems` environment:

```julia
import Pkg
Pkg.activate("/home/lab/NOTEBOOKS/LSHOW/GenLAProblems")
Pkg.develop(path="/home/lab/NOTEBOOKS/LSHOW/LAlatex")
```
