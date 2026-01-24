module GenLAProblems

const _pythoncall_loaded = Ref(false)

_in_precompile() = ccall(:jl_generating_output, Cint, ()) == 1
_pythoncall_disabled() = get(ENV, "GENLAPROBLEMS_DISABLE_PYTHONCALL", "") != ""

function _ensure_pythoncall()
    if _pythoncall_disabled()
        return nothing
    end
    if _in_precompile()
        return nothing
    end
    if !_pythoncall_loaded[]
        @eval using PythonCall
        _pythoncall_loaded[] = true
    end
    return PythonCall
end
using Symbolics
using AbstractAlgebra
import AbstractAlgebra: charpoly
using BlockArrays
using LinearAlgebra
using Random
using Hadamard

using Reexport
@reexport using LAlatex
using PrecompileTools

"""
    Base.adjoint(p::AbstractAlgebra.Generic.Poly{Rational{BigInt}})

Return `p` unchanged to avoid accidental polynomial adjoints.
"""
function Base.adjoint(p::AbstractAlgebra.Generic.Poly{Rational{BigInt}})
    return p
end

"""
    Base.transpose(p::AbstractAlgebra.Generic.Poly{Rational{BigInt}})

Return `p` unchanged to avoid accidental polynomial transposes.
"""
function Base.transpose(p::AbstractAlgebra.Generic.Poly{Rational{BigInt}})
    return p
end

const NO_VALUE = (:none, nothing)

"""
    is_none_val(x) -> Bool

Return `true` when `x` is `:none` or `nothing`.
"""
is_none_val(x) = x === :none || x === nothing

const _la_figures = Ref{Any}(nothing)

struct NMProxy end

const nM = NMProxy()

function (::NMProxy)()
    return nM
end

function _show_svg(svg)
    display(MIME"image/svg+xml"(), svg)
    return svg
end

function Base.getproperty(p::NMProxy, name::Symbol)
    if name === :show_eig_tbl
        return (args...; kwargs...) -> _show_svg(load_la_figures().eig_tbl_svg(args...; kwargs...))
    elseif name === :show_svd_tbl
        return (args...; kwargs...) -> _show_svg(load_la_figures().svd_tbl_svg(args...; kwargs...))
    elseif name === :show_ge_tbl
        return (args...; kwargs...) -> _show_svg(load_la_figures().ge_tbl_svg(args...; kwargs...))
    elseif name === :show_qr_tbl
        return (args...; kwargs...) -> _show_svg(load_la_figures().qr_tbl_svg(args...; kwargs...))
    elseif name === :show_ge
        return (args...; kwargs...) -> _show_svg(load_la_figures().svg(args...; kwargs...))
    elseif name === :show_qr
        return (args...; kwargs...) -> _show_svg(load_la_figures().qr_svg(args...; kwargs...))
    elseif name === :la || name === :la_figures
        return load_la_figures()
    end

    _ensure_pythoncall()
    return getproperty(load_la_figures(), name)
end

"""
    load_la_figures() -> la_figures

Load the Python `la_figures` module via PythonCall.
"""
function load_la_figures()
    if _la_figures[] === nothing
        try
            pc = _ensure_pythoncall()
            if pc === nothing
                return nothing
            end
            _la_figures[] = Base.invokelatest(PythonCall.pyimport, "la_figures")
        catch err
            error(
                "Python module `la_figures` is required by GenLAProblems.\n" *
                "Install it in the active Python environment.\n\n" *
                "Original error:\n$err"
            )
        end
    end
    return _la_figures[]
end


"""
    nM -> NMProxy

Proxy that exposes la_figures display helpers and la_figures attributes.
"""

include("MatrixGeneration.jl")
include("SolveProblems.jl")
include("ge.jl")

export load_la_figures, nM
export symbol_vector, symbols_matrix, form_linear_combination
export invert_unit_lower, unit_lower, lower, gen_full_col_rank_matrix
export ref_matrix, rref_matrix, symmetric_matrix, skew_symmetric_matrix
export e_i, i_with_onecol, gen_permutation_matrix
export W_2_matrix, Q_2_matrix
export W_3_matrix, Q_3_matrix
export Q_4_blocks
export W_4_matrix, Q_4_matrix
export W_matrix, Q_matrix, sparse_W_matrix, sparse_Q_matrix
export split_R_RHS, particular_solution, homogeneous_solutions
export gen_particular_solution
export gen_gj_matrix, gen_rhs, gen_gj_pb
export gen_inv_pb, gen_lu_pb, gen_plu_pb, gen_ldlt_pb
export normal_eq_reduce_to_ref, reduce_to_ref, decorate_ge, ge_variable_type
export ca_projection_matrix
export gen_qr_problem_3, gen_qr_problem_4, gen_qr_problem
export gram_schmidt_w, normalize_columns, qr_layout, gram_schmidt_stable
export gen_eigenproblem, gen_symmetric_eigenproblem, gen_non_diagonalizable_eigenproblem, gen_svd_problem
export gen_cx_eigenproblem
export jordan_block, jordan_form, gen_from_jordan_form, gen_degenerate_matrix
export charpoly
export ge, show_solution
export ShowGe, ref!, show_layout!, show_system, create_cascade!, show_backsubstitution!, show_solution!
export show_backsubstitution, show_forwardsubstitution, solutions
export round_value, round_matrices

# Precompile pure-Julia workloads to reduce latency without PythonCall.
@compile_workload begin
    Random.seed!(1)
    pivot_cols, A = gen_gj_matrix(3, 3, 3)
    gen_rhs(A, pivot_cols)
    ref_matrix(3, 3, 3)
    rref_matrix(3, 3, 3)
    charpoly(A)
    gen_eigenproblem([1, 2, 3])
    gen_symmetric_eigenproblem([1, 2, 3])
    gen_qr_problem(4)
    gen_qr_problem_4()
end

end
