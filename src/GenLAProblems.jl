module GenLAProblems

const __version__ = "1.0.0-DEV"
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
    _ensure_pythoncall()
    builtins = _pyimport("builtins")
    has = Base.invokelatest(PythonCall.pycall, _pygetattr(builtins, :hasattr), obj, String(name))
    if Base.invokelatest(PythonCall.pyconvert, Bool, has)
        return _pygetattr(obj, name)
    end
    sub = _pyimport(mod)
    return _pygetattr(sub, name)
end

"""
    la_version(), la_build(), ml_version(), ml_build()

Return version/build strings exposed by the Python packages `la_figures` and
`matrixlayout`. These require PythonCall at runtime.
"""
function la_version()
    la = load_la_figures()
    v = _pygetattr(la, :__version__)
    return Base.invokelatest(PythonCall.pyconvert, String, v)
end

function la_build()
    la = load_la_figures()
    v = _pygetattr(la, :__build__)
    return Base.invokelatest(PythonCall.pyconvert, String, v)
end

function ml_version()
    ml = load_matrixlayout()
    v = _pygetattr(ml, :__version__)
    return Base.invokelatest(PythonCall.pyconvert, String, v)
end

function ml_build()
    ml = load_matrixlayout()
    v = _pygetattr(ml, :__build__)
    return Base.invokelatest(PythonCall.pyconvert, String, v)
end
import Symbolics
using AbstractAlgebra
import AbstractAlgebra: charpoly
using BlockArrays
using LinearAlgebra
using Random
using Hadamard

using Reexport
@reexport using LAlatex
using LaTeXStrings: LaTeXString
using PrecompileTools
include("SymPyHelpers.jl")
using .SymPyHelpers:
    sym_mat, sym_vec, sym_zero,
    sym_mul, sym_add, sym_pow, sym_eq, sym_is_zero, sym_vec_zero,
    sym_to_julia_vec, sym_to_julia_mat, sym_subs_numeric


"""
    svd_matrices_from_spec(spec; reduced=true)

Build `(U, Σ, V, rank)` from an SVD spec dictionary (as returned by
`nM.show_svd_tbl` / `la_figures.svd_tbl_spec`). The matrix type is preserved:
if the spec contains SymPy objects, returns SymPy matrices; otherwise returns
Julia matrices. When `reduced=true`, only nonzero singular value groups are included.
"""
function svd_matrices_from_spec(spec; reduced::Bool=true)
    _ensure_pythoncall()
    la = load_la_figures()
    svd_from_spec = _pygetattr(la, :svd_matrices_from_spec)
    return _pycall(svd_from_spec, spec; reduced=reduced)
end

"""
    eig_matrices_from_spec(spec; orthonormal=true)

Build `(Λ, V)` from an eigen-table spec dictionary (as returned by
`nM.show_eig_tbl` / `la_figures.eig_tbl_spec`). Uses `qvecs` when available
and `orthonormal=true`, otherwise uses `evecs`. The matrix type is preserved.
"""
function eig_matrices_from_spec(spec; orthonormal::Bool=true)
    _ensure_pythoncall()
    la = load_la_figures()
    eig_from_spec = _pygetattr(la, :eig_matrices_from_spec)
    return _pycall(eig_from_spec, spec; orthonormal=orthonormal)
end

"""
    qr_matrices_from_grid(mats)

Extract QR-related matrices from a Gram–Schmidt QR grid (as returned by
`nM.gram_schmidt_qr`). Returns `(A, W, WtA, WtW, S, Qt, Q, R)` as a named tuple.
The matrix type is preserved (Julia or SymPy).
"""
function qr_matrices_from_grid(mats)
    _ensure_pythoncall()
    la = load_la_figures()
    qr_from_grid = _pygetattr(la, :qr_matrices_from_grid)
    return _pycall(qr_from_grid, mats)
end

"""
    qr_matrices_dict_from_grid(mats)

Return QR matrices as a plain dict for JSON-friendly usage.
"""
function qr_matrices_dict_from_grid(mats)
    _ensure_pythoncall()
    la = load_la_figures()
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

function _show_svg(svg)
    if svg === nothing
        return SVGOut("")
    end
    if _py_is_py(svg)
        _ensure_pythoncall()
        py = getfield(@__MODULE__, :PythonCall)
        if _py_is_none(svg)
            return SVGOut("")
        end
        return SVGOut(Base.invokelatest(py.pyconvert, String, svg))
    elseif svg isa AbstractString
        return SVGOut(svg)
    end
    return svg
end

"""
    py_show_svg(svg)

Display an SVG in a Python notebook (e.g., %%julia cell) via IPython.display.SVG.
Accepts `SVGOut`, raw SVG strings, or PythonCall.Py SVG objects.
"""
function py_show_svg(svg)
    try
        _ensure_pythoncall()
        ip = _pyimport("IPython.display")
        py_display = _pygetattr(ip, :display)
        py_svg = _pygetattr(ip, :SVG)
        if svg isa SVGOut
            return _pycall(py_display, _pycall(py_svg, svg.svg))
        elseif svg isa AbstractString
            return _pycall(py_display, _pycall(py_svg, svg))
        elseif _py_is_py(svg)
            s = Base.invokelatest(getfield(@__MODULE__, :PythonCall).pyconvert, String, svg)
            return _pycall(py_display, _pycall(py_svg, s))
        else
            error("py_show_svg expects SVGOut, SVG string, or PythonCall.Py")
        end
    catch
        # Fallback for Julia kernels without IPython: use Julia's display.
        if svg isa SVGOut
            return Base.display(svg)
        elseif svg isa AbstractString
            return Base.display(SVGOut(svg))
        end
        error("py_show_svg expects SVGOut, SVG string, or PythonCall.Py")
    end
end

"""
    show_svg(svg)

Alias for `py_show_svg`, for notebook-friendly SVG display.
"""
show_svg(svg) = py_show_svg(svg)

"""
    l_show_svd(A, U, Σ, Vt, rankA)

Display an SVD factorization with block structure separating the rank and null
space components.
"""
function l_show_svd(A, U, Σ, Vt, rankA)
    Ub = BlockArray(U, [size(U, 1)], [rankA, size(U, 2) - rankA])
    Σb = BlockArray(Σ, [rankA, size(Σ, 1) - rankA], [rankA, size(Σ, 2) - rankA])
    Vtb = BlockArray(Vt, [rankA, size(Vt, 1) - rankA], [size(Vt, 2)])
    display(l_show(LaTeXString("A = U \\Sigma V^T : "), A, " = ", Ub, Σb, Vtb))
    return nothing
end

function _clean_tmp_kwargs(kwargs)
    clean = Dict(kwargs)
    pop!(clean, :tmp_dir, nothing)
    pop!(clean, :keep_file, nothing)
    pop!(clean, :output_dir, nothing)
    return clean
end

function _map_tmp_to_output(kwargs)
    clean = Dict(kwargs)
    if haskey(clean, :tmp_dir) && !haskey(clean, :output_dir)
        clean[:output_dir] = clean[:tmp_dir]
    end
    pop!(clean, :tmp_dir, nothing)
    pop!(clean, :keep_file, nothing)
    return clean
end

function _normalize_render_opts(render_opts; tmp_dir=nothing)
    opts = if render_opts === nothing
        Dict{String,Any}()
    elseif !(render_opts isa AbstractDict)
        Dict{String,Any}(render_opts)
    else
        Dict{String,Any}(render_opts)
    end
    if tmp_dir !== nothing && !haskey(opts, "output_dir")
        opts["output_dir"] = tmp_dir
    end
    return opts
end

function _bundle_result(dict)
    _ensure_pythoncall()
    py_get = Base.invokelatest(PythonCall.pygetattr, dict, "get")
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
        la = load_la_figures()
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
        la = load_la_figures()
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
    la = load_la_figures()
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
    la = load_la_figures()
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

function _nm_qr_tbl_svg(args...; kwargs...)
    compute_kw, render_kw = _split_qr_tbl_kw(kwargs)
    spec = _qr_tbl_spec_from_args(args...; compute_kw=compute_kw)
    svg = _render_qr_from_spec(spec; render_kw...)
    return svg, spec
end

function Base.getproperty(p::NMProxy, name::Symbol)
    if name === :ge || name === :_to_svg_str
        return matrixlayout_ge
    elseif name === :show_eig_tbl
        return _nm_bundle_wrapper(:eig_tbl_bundle)
    elseif name === :show_svd_tbl
        return _nm_bundle_wrapper(:svd_tbl_bundle)
    elseif name === :show_ge_tbl
        return _nm_svg_wrapper(:ge_tbl_svg; kwcleaner=_clean_tmp_kwargs, fallback_mod="la_figures.ge_convenience")
    elseif name === :show_qr_tbl
        return _nm_bundle_wrapper(:qr_tbl_bundle; kwcleaner=_clean_tmp_kwargs)
    elseif name === :show_ge
        return _nm_svg_wrapper(:svg; kwcleaner=_clean_tmp_kwargs)
    elseif name === :show_qr
        return _nm_svg_wrapper(:qr_svg; kwcleaner=_map_tmp_to_output, with_first_arg=true)
    elseif name === :la || name === :la_figures
        return load_la_figures()
    elseif name === :ml || name === :matrixlayout
        return load_matrixlayout()
    elseif name === :gram_schmidt_qr
        return _nm_gram_schmidt_qr
    elseif name === :qr_tbl_svg
        return _nm_qr_tbl_svg
    elseif name === :qr_svg
        return _nm_svg_wrapper(:qr_svg; kwcleaner=_map_tmp_to_output, wrap_svg=false, with_first_arg=true)
    elseif name === :qr_matrices_from_grid
        return qr_matrices_from_grid
    elseif name === :eig_tbl_svg
        return _nm_bundle_wrapper(:eig_tbl_bundle; wrap_svg=false)
    elseif name === :svd_tbl_svg
        return _nm_bundle_wrapper(:svd_tbl_bundle; wrap_svg=false)
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

export __version__, __build__
export la_version, la_build, ml_version, ml_build
export load_la_figures, load_matrixlayout, nM, sympy
export ensure_pythoncall!
export sym_mat, sym_vec, sym_zero
export sym_mul, sym_add, sym_pow, sym_eq, sym_is_zero, sym_vec_zero
export sym_to_julia_vec, sym_to_julia_mat, sym_subs_numeric
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
export show_ge_final, show_solution, py_show_svg, show_svg, l_show_svd
export ShowGE, ref!, show_layout!, show_system, create_cascade!, show_backsubstitution!, show_solution!
export show_backsubstitution, show_forwardsubstitution, solutions, rhs_block
export round_value, round_matrices
export svd_matrices_from_spec, eig_matrices_from_spec, qr_matrices_from_grid
export qr_matrices_dict_from_grid
export mm_to_px, px_to_mm

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
