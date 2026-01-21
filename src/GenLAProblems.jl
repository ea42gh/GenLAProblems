module GenLAProblems

const _pythoncall_loaded = Ref(false)

function _ensure_pythoncall()
    if !_pythoncall_loaded[]
        @eval using PythonCall
        _pythoncall_loaded[] = true
    end
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

const _itikz = Ref{Any}(nothing)
const _nM = Ref{Any}(nothing)
const _la_figures = Ref{Any}(nothing)
const _nm_proxy = Ref{Any}(nothing)

struct NMProxy
    nm::Any
    la::Any
end

function _show_svg(svg)
    display(MIME"image/svg+xml"(), svg)
    return svg
end

function Base.getproperty(p::NMProxy, name::Symbol)
    if name === :show_eig_tbl
        return (args...; kwargs...) -> _show_svg(p.la.eig_tbl_svg(args...; kwargs...))
    elseif name === :show_svd_tbl
        return (args...; kwargs...) -> _show_svg(p.la.svd_tbl_svg(args...; kwargs...))
    elseif name === :show_ge_tbl
        return (args...; kwargs...) -> _show_svg(p.la.ge_tbl_svg(args...; kwargs...))
    elseif name === :show_qr_tbl
        return (args...; kwargs...) -> _show_svg(p.la.qr_tbl_svg(args...; kwargs...))
    elseif name === :show_ge
        return (args...; kwargs...) -> _show_svg(p.la.svg(args...; kwargs...))
    elseif name === :show_qr
        return (args...; kwargs...) -> _show_svg(p.la.qr_svg(args...; kwargs...))
    elseif name === :la || name === :la_figures
        return p.la
    elseif name === :nm
        return p.nm
    end

    _ensure_pythoncall()
    if p.nm !== nothing && PythonCall.pyhasattr(p.nm, String(name))
        return getproperty(p.nm, name)
    end
    return getproperty(p.la, name)
end

"""
    load_la_figures() -> la_figures

Load the Python `la_figures` module via PythonCall.
"""
function load_la_figures()
    if _la_figures[] === nothing
        try
            _ensure_pythoncall()
            _la_figures[] = pyimport("la_figures")
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
    load_itikz() -> (itikz, nicematrix)

Load the Python `itikz` module and its `nicematrix` helper via PythonCall.
"""
function load_itikz()
    if _itikz[] === nothing
        try
            _ensure_pythoncall()
            _itikz[] = pyimport("itikz")
            _nM[] = pyimport("itikz.nicematrix")
        catch err
            error(
                "Python module `itikz` (and `itikz.nicematrix`) is required by GenLAProblems.\n" *
                "Install it in the active Python environment.\n\n" *
                "Original error:\n$err"
            )
        end
    end
    return _itikz[], _nM[]
end

"""
    nM() -> NMProxy

Return a proxy that exposes la_figures display helpers and forwards other
attributes to `itikz.nicematrix` when available.
"""
function nM()
    if _nm_proxy[] === nothing
        nm = nothing
        try
            _itikz[] = pyimport("itikz")
            nm = pyimport("itikz.nicematrix")
            _nM[] = nm
        catch
            nm = nothing
        end
        _nm_proxy[] = NMProxy(nm, load_la_figures())
    end
    return _nm_proxy[]
end

const nM = nM()

include("MatrixGeneration.jl")
include("SolveProblems.jl")
include("ge.jl")

export load_itikz, load_la_figures, nM
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

end
