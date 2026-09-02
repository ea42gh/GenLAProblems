module GenLAProblems

"""
    __version__::String

Package version string for the loaded `GenLAProblems` module.
"""
const __version__ = "1.1.0"

"""
    __build__::String

Build or source revision metadata for the loaded `GenLAProblems` module.
When unavailable, this is set to `"unknown"`.
"""
const __build__ = "unknown"

using LinearAlgebra
using Random

import Hadamard
using PrecompileTools

const NO_VALUE = (:none, nothing)

"""
    is_none_val(x) -> Bool

Return `true` when `x` is `:none` or `nothing`.
"""
is_none_val(x) = x === :none || x === nothing

include("IntegerFamilies.jl")
include("SymbolNames.jl")
include("BasicMatrices.jl")
include("GEProblems.jl")
include("FactorizationProblems.jl")
include("OrthogonalProblems.jl")
include("SpectralProblems.jl")

export __version__, __build__
export symbol_vector, symbols_matrix
export ca_projection_matrix, gen_matrix
export invert_unit_lower, unit_lower, lower, gen_full_col_rank_matrix
export ref_matrix, rref_matrix, symmetric_matrix, skew_symmetric_matrix
export e_i, i_with_onecol, gen_permutation_matrix
export W_2_matrix, Q_2_matrix
export W_3_matrix, Q_3_matrix
export Q_4_blocks
export W_4_matrix, Q_4_matrix
export W_matrix, Q_matrix, sparse_W_matrix, sparse_Q_matrix
export gen_particular_solution
export gen_gj_matrix, gen_rhs, gen_gj_pb, gen_inconsistent_gj_pb
export gen_inv_pb, gen_lu_pb, gen_plu_pb, gen_ldlt_pb
export gen_qr_problem
export gen_eigenproblem,
    gen_symmetric_eigenproblem, gen_non_diagonalizable_eigenproblem, gen_svd_problem
export gen_cx_eigenproblem
export jordan_block, jordan_form, gen_from_jordan_form, gen_degenerate_matrix

# Precompile pure-Julia workloads to reduce latency without PythonCall.
@compile_workload begin
    Random.seed!(1)
    pivot_cols, A = gen_gj_matrix(3, 3, 3)
    gen_rhs(A, pivot_cols)
    ref_matrix(3, 3, 3)
    rref_matrix(3, 3, 3)
    gen_eigenproblem([1, 2, 3])
    gen_symmetric_eigenproblem([1, 2, 3])
    gen_qr_problem(4; family = :pythagorean)
end

end
