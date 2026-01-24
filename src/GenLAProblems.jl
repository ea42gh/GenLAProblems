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

function _pyimport(name::String)
    _ensure_pythoncall()
    return Base.invokelatest(PythonCall.pyimport, name)
end

function _pycall(f, args...; kwargs...)
    _ensure_pythoncall()
    return Base.invokelatest(PythonCall.pycall, f, args...; kwargs...)
end

function _pygetattr(obj, name::Symbol)
    _ensure_pythoncall()
    return Base.invokelatest(PythonCall.pygetattr, obj, String(name))
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
const _matrixlayout = Ref{Any}(nothing)
const _sympy = Ref{Any}(nothing)

struct SympyProxy end
const sympy = SympyProxy()

struct NMProxy end

const nM = NMProxy()

function (::NMProxy)()
    return nM
end

function _show_svg(svg)
    if _ensure_pythoncall() !== nothing && svg isa PythonCall.Py
        svg = SVGOut(Base.invokelatest(PythonCall.pyconvert, String, svg))
    elseif svg isa AbstractString
        svg = SVGOut(svg)
    end
    return svg
end

function _clean_tmp_kwargs(kwargs)
    clean = Dict(kwargs)
    pop!(clean, :tmp_dir, nothing)
    pop!(clean, :keep_file, nothing)
    pop!(clean, :output_dir, nothing)
    return clean
end

function Base.getproperty(p::NMProxy, name::Symbol)
    if name === :ge || name === :_to_svg_str
        return matrixlayout_ge
    elseif name === :show_eig_tbl
        return function (args...; kwargs...)
            clean = _clean_tmp_kwargs(kwargs)
            la = load_la_figures()
            eig_tbl_svg = _pygetattr(la, :eig_tbl_svg)
            return _show_svg(_pycall(eig_tbl_svg, args...; clean...))
        end
    elseif name === :show_svd_tbl || name === :show_svd_table
        return function (args...; kwargs...)
            clean = _clean_tmp_kwargs(kwargs)
            la = load_la_figures()
            svd_tbl_svg = _pygetattr(la, :svd_tbl_svg)
            return _show_svg(_pycall(svd_tbl_svg, args...; clean...))
        end
    elseif name === :show_ge_tbl
        return function (args...; kwargs...)
            clean = _clean_tmp_kwargs(kwargs)
            la = load_la_figures()
            ge_tbl_svg = _pygetattr(la, :ge_tbl_svg)
            return _show_svg(_pycall(ge_tbl_svg, args...; clean...))
        end
    elseif name === :show_qr_tbl
        return function (args...; kwargs...)
            clean = _clean_tmp_kwargs(kwargs)
            la = load_la_figures()
            qr_tbl_svg = _pygetattr(la, :qr_tbl_svg)
            return _show_svg(_pycall(qr_tbl_svg, args...; clean...))
        end
    elseif name === :show_ge
        return function (args...; kwargs...)
            clean = _clean_tmp_kwargs(kwargs)
            la = load_la_figures()
            svg_fn = _pygetattr(la, :svg)
            return _show_svg(_pycall(svg_fn, args...; clean...))
        end
    elseif name === :show_qr
        return function (args...; kwargs...)
            clean = _clean_tmp_kwargs(kwargs)
            la = load_la_figures()
            qr_svg = _pygetattr(la, :qr_svg)
            return _show_svg(_pycall(qr_svg, args...; clean...))
        end
    elseif name === :la || name === :la_figures
        return load_la_figures()
    elseif name === :ml || name === :matrixlayout
        return load_matrixlayout()
    elseif name === :gram_schmidt_qr
        return function (args...; kwargs...)
            clean = Dict(kwargs)
            if haskey(clean, :tmp_dir) && !haskey(clean, :output_dir)
                clean[:output_dir] = clean[:tmp_dir]
            end
            pop!(clean, :tmp_dir, nothing)
            pop!(clean, :keep_file, nothing)
            la = load_la_figures()
            gram_schmidt_qr = _pygetattr(la, :gram_schmidt_qr)
            svg = _show_svg(_pycall(gram_schmidt_qr, args...; clean...))
            return svg, nothing
        end
    elseif name === :qr_tbl_svg
        return function (args...; kwargs...)
            la = load_la_figures()
            qr_tbl_svg = _pygetattr(la, :qr_tbl_svg)
            return _pycall(qr_tbl_svg, args...; kwargs...)
        end
    elseif name === :qr_svg
        return function (args...; kwargs...)
            la = load_la_figures()
            qr_svg = _pygetattr(la, :qr_svg)
            return _pycall(qr_svg, args...; kwargs...)
        end
    elseif name === :eig_tbl_svg
        return function (args...; kwargs...)
            la = load_la_figures()
            eig_tbl_svg = _pygetattr(la, :eig_tbl_svg)
            return _pycall(eig_tbl_svg, args...; kwargs...)
        end
    elseif name === :svd_tbl_svg
        return function (args...; kwargs...)
            la = load_la_figures()
            svd_tbl_svg = _pygetattr(la, :svd_tbl_svg)
            return _pycall(svd_tbl_svg, args...; kwargs...)
        end
    end

    _ensure_pythoncall()
    return getproperty(load_matrixlayout(), name)
end

function Base.getproperty(::SympyProxy, name::Symbol)
    if _sympy[] === nothing
        _sympy[] = _pyimport("sympy")
    end
    attr = getproperty(_sympy[], name)
    builtins = _pyimport("builtins")
    if PythonCall.pyconvert(Bool, _pycall(builtins.callable, attr))
        return (args...; kwargs...) -> _pycall(attr, args...; kwargs...)
    end
    return attr
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
    load_matrixlayout() -> matrixlayout

Load the Python `matrixlayout` module via PythonCall.
"""
function load_matrixlayout()
    if _matrixlayout[] === nothing
        try
            pc = _ensure_pythoncall()
            if pc === nothing
                return nothing
            end
            _matrixlayout[] = Base.invokelatest(PythonCall.pyimport, "matrixlayout")
        catch err
            error(
                "Python module `matrixlayout` is required by GenLAProblems.\n" *
                "Install it in the active Python environment.\n\n" *
                "Original error:\n$err"
            )
        end
    end
    return _matrixlayout[]
end


"""
    nM -> NMProxy

Proxy that exposes matrixlayout helpers by default and la_figures display helpers.
"""

include("MatrixGeneration.jl")
include("SolveProblems.jl")
include("ge.jl")

export load_la_figures, load_matrixlayout, nM, sympy
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
