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

For the full exported API reference and generator notes, see the
[Usage and Quickstart](https://ea42gh.github.io/GenLAProblems.jl/dev/usage/) page.

## Installation

```julia
import Pkg
Pkg.add("GenLAProblems")
```

For notebook-friendly matrix display, install `LAlatex` separately:

```julia
import Pkg
Pkg.add("LAlatex")
```

## Minimal example

```julia
using GenLAProblems

A, X, B = gen_gj_pb(3, 4, 3; maxint=3, num_rhs=1)
A_bad, B_bad = gen_inconsistent_gj_pb(4, 6, 3; maxint=2, num_rhs=1)
```

With `LAlatex` installed, the same example can be rendered in a notebook:

```julia
using GenLAProblems, LAlatex

A, X, B = gen_gj_pb(3, 4, 3; maxint=3, num_rhs=1)
l_show("A X = B : ", A, X, " = ", B)
```

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

## Intended audience

GenLAProblems is for:

- instructors building worksheets, notebooks, and demos
- students exploring exact linear algebra examples
- developers who need small exact matrices with controlled algebraic structure

It is not intended to replace general-purpose numerical linear algebra libraries.

## Generator notes

- `gen_gj_pb(...)` and `gen_inconsistent_gj_pb(...)` are the canonical exact
  generators for consistent and inconsistent Gaussian-elimination /
  Gauss-Jordan system problems.
- `W_matrix` is the general exact `W` constructor, and `W_2_matrix`,
  `W_3_matrix`, and `W_4_matrix` are the specialized small-size families.
  All consistently return `(d, A)` where `A' * A` is diagonal.
- `Q_matrix` is the general exact orthogonal constructor, and `Q_2_matrix`,
  `Q_3_matrix`, and `Q_4_matrix` are the specialized small-size families.
- `gen_qr_problem(n; family=:auto, maxint=3)` is now the canonical QR-problem
  entrypoint. Supported families are `:pythagorean`, `:hadamard`, `:cayley`,
  and `:sparse`.
- `gen_eigenproblem` and `gen_symmetric_eigenproblem` accept either a standard
  vector of eigenvalues or a `1×n` / `n×1` matrix input.
- `gen_svd_problem(m, n, σ; left_family=:sparse, right_family=:sparse, maxint=3)`
  now lets you choose the orthogonal-family construction for the left and right
  SVD factors independently.
  The `:hadamard` family is exact-rational only when the size both supports a
  Hadamard matrix and has an integer square root (for example `n = 4`, but not
  `n = 12`).

## Package boundaries

GenLAProblems is intentionally generator-focused.

- Use `GenLAProblems` to construct exact matrices and exercise instances.
- Import `LAlatex` separately when you want notebook-friendly LaTeX display.
- Import `LATeachingSuite` when you want GE reduction workflows, `ShowGE`,
  back-substitution displays, or Python-backed figure helpers.
- The Julia-facing Python bridge is owned by `LATeachingSuite.PythonBridge`,
  not by `GenLAProblems`.
- Removed render/display APIs such as the old `nM` proxy and bundle/spec
  extraction helpers are no longer part of `GenLAProblems`; use
  `LATeachingSuite` for the current workflow and rendering helpers.

The package itself no longer depends on `LAlatex`; display support is an
explicit companion-layer choice rather than a hard dependency of the generator
core.

## Local development

If you are working from this monorepo, install the local `LAlatex` checkout in the
`GenLAProblems` environment:

```julia
import Pkg
Pkg.activate("/home/lab/NOTEBOOKS/LA/GenLAProblems")
Pkg.develop(path="/home/lab/NOTEBOOKS/LA/LAlatex")
```
