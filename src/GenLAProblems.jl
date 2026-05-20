module GenLAProblems

const __version__ = "1.0.0"
const __build__ = "908575d 2026-02-06T15:37:28-05:00"

const _pythoncall_loaded = Ref(false)
const _pythoncall_mod = Ref{Any}(nothing)

_in_precompile() = ccall(:jl_generating_output, Cint, ()) == 1
_pythoncall_disabled() = get(ENV, "GENLAPROBLEMS_DISABLE_PYTHONCALL", "") != ""

function _ensure_pythoncall()
    if _pythoncall_disabled()
        throw(ArgumentError("PythonCall disabled. Unset GENLAPROBLEMS_DISABLE_PYTHONCALL or set JULIA_PYTHONCALL_EXE."))
    end
    if _in_precompile()
        throw(ArgumentError("PythonCall unavailable during precompile. Call from runtime context."))
    end
    if !_pythoncall_loaded[]
        @eval using PythonCall
        _pythoncall_loaded[] = true
    end
    if _pythoncall_mod[] === nothing
        _pythoncall_mod[] = Base.invokelatest(() -> PythonCall)
    end
    return _pythoncall_mod[]
end

"""
    ensure_pythoncall!()

Explicitly initialize PythonCall and return the module handle. Raises a
clear error when PythonCall is disabled or unavailable.
"""
ensure_pythoncall!() = _ensure_pythoncall()

function _pyimport(name::String)
    py = _ensure_pythoncall()
    return Base.invokelatest(py.pyimport, name)
end

function _pycall(f, args...; kwargs...)
    py = _ensure_pythoncall()
    return Base.invokelatest(py.pycall, f, args...; kwargs...)
end

function _pygetattr(obj, name::Symbol)
    py = _ensure_pythoncall()
    return Base.invokelatest(py.pygetattr, obj, String(name))
end

function _pygetattr_fallback(obj, name::Symbol, mod::String)
    py = _ensure_pythoncall()
    builtins = _pyimport("builtins")
    has = Base.invokelatest(py.pycall, _pygetattr(builtins, :hasattr), obj, String(name))
    if Base.invokelatest(py.pyconvert, Bool, has)
        return _pygetattr(obj, name)
    end
    sub = _pyimport(mod)
    return _pygetattr(sub, name)
end

using LinearAlgebra
using Random
using AbstractAlgebra
import AbstractAlgebra: charpoly

using Reexport
using LAlatex
using LaTeXStrings: LaTeXString
using PrecompileTools
include("SymPyHelpers.jl")
using .SymPyHelpers:
    sym_mat, sym_vec, sym_zero,
    sym_mul, sym_add, sym_pow, sym_eq, sym_is_zero, sym_vec_zero,
    sym_to_julia_vec, sym_to_julia_mat, sym_subs_numeric

const _symbolics_mod = Ref{Any}(nothing)
const _blockarrays_mod = Ref{Any}(nothing)
const _hadamard_mod = Ref{Any}(nothing)

function _ensure_symbolics()
    if _symbolics_mod[] === nothing
        @eval import Symbolics
        _symbolics_mod[] = Base.invokelatest(() -> Symbolics)
    end
    return _symbolics_mod[]
end

function _ensure_blockarrays()
    if _blockarrays_mod[] === nothing
        @eval import BlockArrays
        _blockarrays_mod[] = Base.invokelatest(() -> BlockArrays)
    end
    return _blockarrays_mod[]
end

function _ensure_hadamard()
    if _hadamard_mod[] === nothing
        @eval import Hadamard
        _hadamard_mod[] = Base.invokelatest(() -> Hadamard)
    end
    return _hadamard_mod[]
end


"""
    svd_matrices_from_spec(spec; reduced=true)

Build `(U, Σ, V, rank)` from an SVD spec dictionary (as returned by
`LATeachingSuite.svd_bundle` / `LAFigureSpecs.svd_tbl_spec`). The matrix type is preserved:
if the spec contains SymPy objects, returns SymPy matrices; otherwise returns
Julia matrices. When `reduced=true`, only nonzero singular value groups are included.
"""
function svd_matrices_from_spec(spec; reduced::Bool=true)
    _ensure_pythoncall()
    la = load_LAFigureSpecs()
    svd_from_spec = _pygetattr(la, :svd_matrices_from_spec)
    return materialize_python_value(_pycall(svd_from_spec, spec; reduced=reduced))
end

"""
    eig_matrices_from_spec(spec; orthonormal=true)

Build `(Λ, V)` from an eigen-table spec dictionary (as returned by
`LATeachingSuite.eig_bundle` / `LAFigureSpecs.eig_tbl_spec`). Uses `qvecs` when available
and `orthonormal=true`, otherwise uses `evecs`. The matrix type is preserved.
"""
function eig_matrices_from_spec(spec; orthonormal::Bool=true)
    _ensure_pythoncall()
    la = load_LAFigureSpecs()
    eig_from_spec = _pygetattr(la, :eig_matrices_from_spec)
    return materialize_python_value(_pycall(eig_from_spec, spec; orthonormal=orthonormal))
end

"""
    qr_matrices_from_grid(mats)

Extract QR-related matrices from a Gram–Schmidt QR grid (as returned by
`LATeachingSuite.qr_figure`). Returns `(A, W, WtA, WtW, S, Qt, Q, R)` as a named tuple.
The matrix type is preserved (Julia or SymPy).
"""
function qr_matrices_from_grid(mats)
    _ensure_pythoncall()
    la = load_LAFigureSpecs()
    qr_from_grid = _pygetattr(la, :qr_matrices_from_grid)
    qr = _pycall(qr_from_grid, mats)
    getmat(name::String) = begin
        try
            materialize_python_value(qr[name])
        catch
            nothing
        end
    end
    return (
        A = getmat("A"),
        W = getmat("W"),
        WtA = getmat("WtA"),
        WtW = getmat("WtW"),
        S = getmat("S"),
        Qt = getmat("Qt"),
        Q = getmat("Q"),
        R = getmat("R"),
    )
end

"""
    qr_matrices_dict_from_grid(mats)

Return QR matrices as a plain dict for JSON-friendly usage.
"""
function qr_matrices_dict_from_grid(mats)
    _ensure_pythoncall()
    la = load_LAFigureSpecs()
    qr_from_grid = _pygetattr(la, :qr_matrices_dict_from_grid)
    return _pycall(qr_from_grid, mats)
end

const NO_VALUE = (:none, nothing)

"""
    is_none_val(x) -> Bool

Return `true` when `x` is `:none` or `nothing`.
"""
is_none_val(x) = x === :none || x === nothing

"""
    mm_to_px(mm) -> Float64

Convert millimeters to px-equivalent SVG units (96 px per inch).
"""
mm_to_px(mm::Real) = float(mm) * 96.0 / 25.4

"""
    px_to_mm(px) -> Float64

Convert px-equivalent SVG units to millimeters (96 px per inch).
"""
px_to_mm(px::Real) = float(px) * 25.4 / 96.0

const _LAFigureSpecs = Ref{Any}(nothing)
const _matrixlayout = Ref{Any}(nothing)
const _sympy = Ref{Any}(nothing)

function _py_is_none(svg)
    if !isdefined(@__MODULE__, :PythonCall)
        return false
    end
    try
        py = getfield(@__MODULE__, :PythonCall)
        return Base.invokelatest(py.pyconvert, Any, svg) === nothing
    catch
        return false
    end
end

function _py_is_py(svg)
    if !isdefined(@__MODULE__, :PythonCall)
        return false
    end
    py = getfield(@__MODULE__, :PythonCall)
    return svg isa py.Py
end

include("PythonBridgeUtils.jl")
include("DisplayInterop.jl")

function _bundle_result(dict)
    py = _ensure_pythoncall()
    py_get = Base.invokelatest(py.pygetattr, dict, "get")
    spec = _pycall(py_get, "spec")
    svg = _pycall(py_get, "svg")
    if _py_is_py(svg) && _py_is_none(svg)
        svg = nothing
    end
    return spec, svg
end

function _nm_bundle_wrapper(bundle_sym::Symbol; kwcleaner=_map_tmp_to_output, wrap_svg::Bool=true)
    return function (args...; kwargs...)
        clean = kwcleaner(kwargs)
        la = load_LAFigureSpecs()
        bundle_fn = _pygetattr(la, bundle_sym)
        spec, svg = _bundle_result(_pycall(bundle_fn, args...; clean...))
        return wrap_svg ? (_show_svg(svg), spec) : (svg, spec)
    end
end

function _nm_svg_wrapper(svg_sym::Symbol;
    kwcleaner=_clean_tmp_kwargs,
    fallback_mod::Union{Nothing,String}=nothing,
    wrap_svg::Bool=true,
    with_first_arg::Bool=false,
)
    return function (args...; kwargs...)
        clean = kwcleaner(kwargs)
        la = load_LAFigureSpecs()
        svg_fn = fallback_mod === nothing ? _pygetattr(la, svg_sym) : _pygetattr_fallback(la, svg_sym, fallback_mod)
        svg = _pycall(svg_fn, args...; clean...)
        if wrap_svg
            return with_first_arg ? (_show_svg(svg), args[1]) : _show_svg(svg)
        end
        return with_first_arg ? (svg, args[1]) : svg
    end
end

function _split_qr_kw(kwargs)
    # Separate kwargs for:
    # - gram_schmidt_qr_matrices (compute)
    # - qr_tbl_spec_from_matrices (spec)
    # - render_qr_svg (render)
    matrices_keys = Set([:allow_rank_deficient, :rank_deficient])
    spec_keys = Set([
        :array_names,
        :fig_scale,
        :preamble,
        :extension,
        :nice_options,
        :label_color,
        :label_text_color,
        :known_zero_color,
        :decorators,
        :strict,
    ])
    render_keys = Set([
        :formatter,
        :toolchain_name,
        :crop,
        :padding,
        :frame,
        :exact_bbox,
        :output_dir,
        :output_stem,
        :tmp_dir,
        :render_opts,
        :strict,
    ])
    matrices_kw = Dict()
    spec_kw = Dict()
    render_kw = Dict()
    for (k, v) in kwargs
        if k === :strict
            spec_kw[k] = v
            render_kw[k] = v
        elseif k in render_keys
            render_kw[k] = v
        elseif k in spec_keys
            spec_kw[k] = v
        elseif k in matrices_keys
            matrices_kw[k] = v
        else
            spec_kw[k] = v
        end
    end
    if haskey(render_kw, :tmp_dir) && !haskey(render_kw, :output_dir)
        render_kw[:output_dir] = render_kw[:tmp_dir]
    end
    pop!(render_kw, :tmp_dir, nothing)
    pop!(render_kw, :keep_file, nothing)
    return matrices_kw, spec_kw, render_kw
end

function _split_qr_tbl_kw(kwargs)
    compute_keys = Set([
        :array_names,
        :fig_scale,
        :preamble,
        :extension,
        :nice_options,
        :label_color,
        :label_text_color,
        :known_zero_color,
        :decorators,
        :strict,
    ])
    render_keys = Set([
        :formatter,
        :toolchain_name,
        :crop,
        :padding,
        :frame,
        :exact_bbox,
        :output_dir,
        :output_stem,
        :tmp_dir,
        :render_opts,
        :strict,
    ])
    compute_kw = Dict()
    render_kw = Dict()
    for (k, v) in kwargs
        if k === :strict
            compute_kw[k] = v
            render_kw[k] = v
        elseif k in render_keys
            render_kw[k] = v
        elseif k in compute_keys
            compute_kw[k] = v
        else
            compute_kw[k] = v
        end
    end
    if haskey(render_kw, :tmp_dir) && !haskey(render_kw, :output_dir)
        render_kw[:output_dir] = render_kw[:tmp_dir]
    end
    pop!(render_kw, :tmp_dir, nothing)
    pop!(render_kw, :keep_file, nothing)
    return compute_kw, render_kw
end

function _qr_spec_from_args(args...; matrices_kw=Dict(), spec_kw=Dict())
    if length(args) != 1
        throw(ArgumentError("qr routines expect a single matrix A; W is computed internally"))
    end
    la = load_LAFigureSpecs()
    gram_schmidt_qr_matrices = _pygetattr(la, :gram_schmidt_qr_matrices)
    qr_tbl_spec_from_matrices = _pygetattr(la, :qr_tbl_spec_from_matrices)
    matrices_nt = (; matrices_kw...)
    spec_nt = (; spec_kw...)
    matrices = _pycall(gram_schmidt_qr_matrices, args...; matrices_nt...)
    spec = _pycall(qr_tbl_spec_from_matrices, matrices; spec_nt...)
    return spec, matrices
end

function _qr_tbl_spec_from_args(args...; compute_kw=Dict())
    if length(args) != 1
        throw(ArgumentError("qr routines expect a single matrix A; W is computed internally"))
    end
    la = load_LAFigureSpecs()
    qr_tbl_spec = _pygetattr(la, :qr_tbl_spec)
    compute_nt = (; compute_kw...)
    return _pycall(qr_tbl_spec, args...; compute_nt...)
end

function _render_qr_from_spec(spec; render_kw...)
    render_qr_svg = _pygetattr(load_matrixlayout(), :render_qr_svg)
    return _pycall(render_qr_svg; spec=spec, render_kw...)
end

function _nm_gram_schmidt_qr(args...; kwargs...)
    matrices_kw, spec_kw, render_kw = _split_qr_kw(kwargs)
    spec, matrices = _qr_spec_from_args(args...; matrices_kw=matrices_kw, spec_kw=spec_kw)
    svg = _render_qr_from_spec(spec; render_kw...)
    return _show_svg(svg), matrices
end

function _nm_qr_svg(args...; kwargs...)
    return _nm_svg_wrapper(:qr_svg; kwcleaner=_map_tmp_to_output)(args...; kwargs...)
end

function _nm_ge_svg(args...; kwargs...)
    return matrixlayout_ge(args...; kwargs...)
end

function _nm_qr_tbl_svg(args...; kwargs...)
    compute_kw, render_kw = _split_qr_tbl_kw(kwargs)
    spec = _qr_tbl_spec_from_args(args...; compute_kw=compute_kw)
    svg = _render_qr_from_spec(spec; render_kw...)
    return svg, spec
end

function _nm_depwarn(name::Symbol, replacement::AbstractString)
    Base.depwarn("`nM.$name` is deprecated; use `$replacement` instead.", Symbol("nM_$(name)"))
end

function Base.getproperty(p::NMProxy, name::Symbol)
    if name === :ge || name === :_to_svg_str
        _nm_depwarn(name, "LATeachingSuite.ge_svg")
        return _nm_ge_svg
    elseif name === :show_eig_tbl
        _nm_depwarn(name, "LATeachingSuite.eig_bundle")
        return _nm_bundle_wrapper(:eig_bundle)
    elseif name === :show_svd_tbl
        _nm_depwarn(name, "LATeachingSuite.svd_bundle")
        return _nm_bundle_wrapper(:svd_bundle)
    elseif name === :show_ge_tbl
        _nm_depwarn(name, "LATeachingSuite.ge_svg (or LATeachingSuite.ge_bundle if you also need the spec)")
        return _nm_svg_wrapper(:ge_tbl_svg; kwcleaner=_clean_tmp_kwargs, fallback_mod="LAFigureSpecs.ge_convenience")
    elseif name === :show_qr_tbl
        _nm_depwarn(name, "LATeachingSuite.qr_bundle")
        return _nm_bundle_wrapper(:qr_bundle; kwcleaner=_clean_tmp_kwargs)
    elseif name === :show_ge
        _nm_depwarn(name, "LATeachingSuite.ge_svg")
        return _nm_svg_wrapper(:ge_svg; kwcleaner=_clean_tmp_kwargs)
    elseif name === :show_qr
        _nm_depwarn(name, "LATeachingSuite.qr_svg (or LATeachingSuite.qr_figure if you also need the computed matrices)")
        return _nm_svg_wrapper(:qr_svg; kwcleaner=_map_tmp_to_output, with_first_arg=true)
    elseif name === :la || name === :LAFigureSpecs
        _nm_depwarn(name, "LATeachingSuite.load_LAFigureSpecs()")
        return load_LAFigureSpecs()
    elseif name === :ml || name === :matrixlayout
        _nm_depwarn(name, "LATeachingSuite.load_matrixlayout()")
        return load_matrixlayout()
    elseif name === :gram_schmidt_qr
        _nm_depwarn(name, "LATeachingSuite.qr_figure")
        return _nm_gram_schmidt_qr
    elseif name === :qr_tbl_svg
        _nm_depwarn(name, "LATeachingSuite.qr_svg")
        return _nm_qr_tbl_svg
    elseif name === :qr_svg
        _nm_depwarn(name, "LATeachingSuite.qr_svg")
        return _nm_svg_wrapper(:qr_svg; kwcleaner=_map_tmp_to_output, wrap_svg=false, with_first_arg=true)
    elseif name === :qr_matrices_from_grid
        _nm_depwarn(name, "LATeachingSuite.qr_matrices_from_grid")
        return qr_matrices_from_grid
    elseif name === :eig_tbl_svg
        _nm_depwarn(name, "LATeachingSuite.eig_svg")
        return _nm_bundle_wrapper(:eig_bundle; wrap_svg=false)
    elseif name === :svd_tbl_svg
        _nm_depwarn(name, "LATeachingSuite.svd_svg")
        return _nm_bundle_wrapper(:svd_bundle; wrap_svg=false)
    end

    _ensure_pythoncall()
    return getproperty(load_matrixlayout(), name)
end

function Base.getproperty(::SympyProxy, name::Symbol)
    if _sympy[] === nothing
        _sympy[] = _pyimport("sympy")
    end
    attr = _pygetattr(_sympy[], name)
    builtins = _pyimport("builtins")
    py = _ensure_pythoncall()
    if Base.invokelatest(py.pyconvert, Bool, _pycall(builtins.callable, attr))
        return (args...; kwargs...) -> _pycall(attr, args...; kwargs...)
    end
    return attr
end

include("MatrixGeneration.jl")
include("SolveProblems.jl")

module PythonBridge

using Reexport

@reexport using ..GenLAProblems:
    ensure_pythoncall!,
    load_LAFigureSpecs,
    load_matrixlayout,
    la_version,
    la_build,
    ml_version,
    ml_build,
    materialize_python_value,
    sympy,
    svd_matrices_from_spec,
    eig_matrices_from_spec,
    qr_matrices_from_grid,
    qr_matrices_dict_from_grid

end

export __version__, __build__
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
export gen_eigenproblem, gen_symmetric_eigenproblem, gen_non_diagonalizable_eigenproblem, gen_svd_problem
export gen_cx_eigenproblem
export jordan_block, jordan_form, gen_from_jordan_form, gen_degenerate_matrix

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
    gen_qr_problem(4; family=:pythagorean)
end

end
