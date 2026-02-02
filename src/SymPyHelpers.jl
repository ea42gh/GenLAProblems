module SymPyHelpers

export sym_mat, sym_vec, sym_zero, sym_mul, sym_add, sym_pow, sym_eq, sym_is_zero, sym_vec_zero
export sym_to_julia_vec, sym_to_julia_mat, sym_subs_numeric

using PythonCall
using ..GenLAProblems: import_sympy

const _sympy = Ref{Any}(nothing)

function _sympy_module()
    if _sympy[] === nothing
        _sympy[] = import_sympy()
    end
    return _sympy[]
end

sym_mat(x) = x isa PythonCall.Py ? x : _sympy_module().Matrix(x)
sym_vec(x) = x isa PythonCall.Py ? x : _sympy_module().Matrix(x)
sym_zero() = _sympy_module().Integer(0)

sym_mul(A, v) = sym_mat(A) * sym_vec(v)
sym_add(A, B) = sym_mat(A) + sym_mat(B)
sym_pow(A, k) = sym_mat(A) ^ k

sym_is_zero(x) = PythonCall.pyconvert(Bool, _sympy_module().simplify(x).is_zero_matrix)
sym_eq(A, B) = sym_is_zero(sym_mat(A) - sym_mat(B))

sym_vec_zero(v) = all(PythonCall.pyconvert(Bool, _sympy_module().simplify(e) == 0) for e in v)

sym_to_julia_vec(x) = x isa PythonCall.Py ? Base.invokelatest(PythonCall.pyconvert, Vector{Any}, x) : x
sym_to_julia_mat(x) = x isa PythonCall.Py ? Base.invokelatest(PythonCall.pyconvert, Matrix{Any}, x) : x

function _sympy_scalar_to_julia(x)
    if !(x isa PythonCall.Py)
        return x
    end
    sympy = _sympy_module()
    is_int = Base.invokelatest(PythonCall.pyconvert, Bool, Base.invokelatest(PythonCall.pygetattr, x, "is_Integer"))
    if is_int
        return Base.invokelatest(PythonCall.pyconvert, Int, x)
    end
    is_rat = Base.invokelatest(PythonCall.pyconvert, Bool, Base.invokelatest(PythonCall.pygetattr, x, "is_Rational"))
    if is_rat
        p = Base.invokelatest(PythonCall.pyconvert, Int, Base.invokelatest(PythonCall.pygetattr, x, "p"))
        q = Base.invokelatest(PythonCall.pyconvert, Int, Base.invokelatest(PythonCall.pygetattr, x, "q"))
        return p//q
    end
    return Base.invokelatest(PythonCall.pyconvert, Float64, sympy.N(x))
end

function _promote_matrix(M::AbstractArray)
    types = Set{DataType}()
    for v in M
        push!(types, typeof(v))
    end
    T = foldl(promote_type, collect(types); init=Any)
    return Array{T}(M)
end

"""
    sym_subs_numeric(A, subs) -> Union{Py, AbstractArray}

Substitute `subs` into a SymPy matrix or Julia matrix convertible to SymPy.
`subs` can be a Dict or list of pairs. Returns a SymPy matrix if symbols remain;
otherwise returns a Julia numeric array.
"""
function sym_subs_numeric(A, subs)
    sympy = _sympy_module()
    symA = sym_mat(A)
    sub_list = subs isa AbstractDict ? collect(pairs(subs)) : subs
    if sub_list isa Pair
        sub_list = [sub_list]
    end
    if sub_list isa AbstractVector
        sub_list = [(p.first, p.second) for p in sub_list]
    end
    subbed = symA.subs(sub_list)
    free = Base.invokelatest(PythonCall.pygetattr, subbed, "free_symbols")
    blen = Base.invokelatest(PythonCall.pybuiltins.len, free)
    nfree = Base.invokelatest(PythonCall.pyconvert, Int, blen)
    if nfree != 0
        return subbed
    end
    M = sym_to_julia_mat(subbed)
    try
        num = map(_sympy_scalar_to_julia, M)
        return _promote_matrix(num)
    catch
        return subbed
    end
end

# Backward-compatible aliases (prefer sym_* names).
const sympy_mat = sym_mat
const sympy_vec = sym_vec
const sympy_zero = sym_zero
const sympy_to_julia_vec = sym_to_julia_vec
const sympy_to_julia_mat = sym_to_julia_mat
const sympy_subs_numeric = sym_subs_numeric

end
