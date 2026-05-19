#!/usr/bin/env bash
set -euo pipefail

cd /home/jovyan/GenLAProblems

report=/home/jovyan/work/kernel_probe_report.txt
: > "$report"

run_probe() {
  local label="$1"
  local code="$2"
  local start end status
  {
    echo "=== $label ==="
    start=$(date +%s)
    set +e
    julia --project=. -e "$code"
    status=$?
    set -e
    end=$(date +%s)
    echo "elapsed=$((end-start))s exit=$status"
    echo
  } >>"$report" 2>&1
}

run_probe "using Random" 'using Random; Random.seed!(42); println("ok")'
run_probe "using LinearAlgebra" 'using LinearAlgebra; println("ok")'
run_probe "using LaTeXStrings" 'using LaTeXStrings; println("ok")'
run_probe "using LAlatex" 'using LAlatex; println("ok")'
run_probe "using GenLAProblems" 'using GenLAProblems; println("ok")'
run_probe "notebook first cell" 'using GenLAProblems, LAlatex, LinearAlgebra, LaTeXStrings, Random; Random.seed!(42); println("ok")'
run_probe "notebook warmup script" 'include("binder/notebook_warmup.jl"); println("ok")'
run_probe "gen_gj_pb" 'using GenLAProblems; A,X,B = gen_gj_pb(3,4,3; maxint=3, num_rhs=2); println(size(A), size(X), size(B))'
run_probe "LAlatex display" 'using GenLAProblems, LAlatex, LaTeXStrings; A,X,B = gen_gj_pb(3,4,3; maxint=3, num_rhs=2); display(LAlatex.l_show(L"A = ", A)); println("ok")'
