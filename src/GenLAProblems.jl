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
using PrecompileTools
include("SymPyHelpers.jl")
using .SymPyHelpers:
    sym_mat, sym_vec, sym_zero,
    sym_mul, sym_add, sym_pow, sym_eq, sym_is_zero, sym_vec_zero,
    sym_to_julia_vec, sym_to_julia_mat, sym_subs_numeric

const _symbolics_mod = Ref{Any}(nothing)
const _hadamard_mod = Ref{Any}(nothing)

function _ensure_symbolics()
    if _symbolics_mod[] === nothing
        @eval import Symbolics
        _symbolics_mod[] = Base.invokelatest(() -> Symbolics)
    end
    return _symbolics_mod[]
end

function _ensure_hadamard()
    if _hadamard_mod[] === nothing
        @eval import Hadamard
        _hadamard_mod[] = Base.invokelatest(() -> Hadamard)
    end
    return _hadamard_mod[]
end


const NO_VALUE = (:none, nothing)

"""
    is_none_val(x) -> Bool

Return `true` when `x` is `:none` or `nothing`.
"""
is_none_val(x) = x === :none || x === nothing

const _sympy = Ref{Any}(nothing)

include("PythonBridgeUtils.jl")

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
