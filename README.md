# GenLAProblems

[![CI](https://github.com/ea42gh/GenLAProblems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ea42gh/GenLAProblems.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ea42gh/GenLAProblems.jl/HEAD?filepath=GenProblems.ipynb)
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://ea42gh.github.io/GenLAProblems.jl/)

## Overview

GenLAProblems is the exercise-generation layer of this teaching stack.
It provides exact linear algebra problem generators and supporting matrix
constructors for standard classroom topics such as GE/GJ, LU/PLU/LDLT, QR,
eigenvalue problems, and SVD.

It is not a general-purpose numerical linear algebra package. The focus is on
exact, classroom-friendly examples rather than large-scale numerical
computation.

The Binder notebook opens `GenProblems.ipynb` and shows the current
generator-focused workflow.

## Why This Package Exists

GenLAProblems exists to make it easy to generate linear algebra exercises whose
arithmetic is still readable by students and instructors.

- It favors exact arithmetic over floating-point approximations.
- It generates standard exercise families used in linear algebra courses.
- It is designed for notebooks, worksheets, and teaching demos rather than
  benchmarking or production numerical workflows.
- It pairs naturally with:
  - `LAlatex` for matrix/factorization display
  - `LATeachingSuite` for Julia-side reduction workflows and display helpers
  - `LAFigureSpecs` and `matrixlayout` for the Python-backed figure/rendering stack

## Generator updates

Recent generator changes to note:

- `gen_qr_problem(n; family=:auto, maxint=3)` is now the canonical QR-problem
  entrypoint. Supported families are `:pythagorean`, `:hadamard`, `:cayley`,
  and `:sparse`.
- `gen_svd_problem(m, n, σ; left_family=:sparse, right_family=:sparse, maxint=3)`
  now lets you choose the orthogonal-family construction for the left and right
  SVD factors independently.
  The `:hadamard` family is exact-rational only when the size both supports a
  Hadamard matrix and has an integer square root (for example `n = 4`, but not
  `n = 12`).
- `W_2_matrix`, `W_3_matrix`, `W_4_matrix`, and `W_matrix` consistently return
  `(d, A)` where `A' * A` is diagonal. `Q_2_matrix`, `Q_3_matrix`, and
  `Q_4_matrix` return orthogonal matrices.
- `gen_eigenproblem` and `gen_symmetric_eigenproblem` accept either a standard
  vector of eigenvalues or a `1×n` / `n×1` matrix input.

## Package boundaries

GenLAProblems is intentionally generator-focused.

- Use `GenLAProblems` to construct exact matrices and exercise instances.
- Import `LAlatex` separately when you want notebook-friendly LaTeX display
  (`Pkg.add("LAlatex")` works directly from the General registry).
- Import `LATeachingSuite` when you want GE reduction workflows, `ShowGE`,
  back-substitution displays, or Python-backed figure helpers.

## Local development

If you are working from this monorepo, install the local `LAlatex` checkout in the
`GenLAProblems` environment:

```julia
import Pkg
Pkg.activate("/home/lab/NOTEBOOKS/LA/GenLAProblems")
Pkg.develop(path="/home/lab/NOTEBOOKS/LA/LAlatex")
```
