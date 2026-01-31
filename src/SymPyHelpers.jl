module SymPyHelpers

export sympy_mat, sympy_vec, sympy_zero, sym_mul, sym_add, sym_pow, sym_eq, sym_is_zero, sym_vec_zero
export sympy_to_julia_vec, sympy_to_julia_mat, sympy_subs_numeric

using PythonCall
using ..GenLAProblems: import_sympy

const _sympy = Ref{Any}(nothing)

function _sympy_module()
    if _sympy[] === nothing
        _sympy[] = import_sympy()
    end
    return _sympy[]
end

sympy_mat(x) = x isa PythonCall.Py ? x : _sympy_module().Matrix(x)
sympy_vec(x) = x isa PythonCall.Py ? x : _sympy_module().Matrix(x)
sympy_zero() = _sympy_module().Integer(0)

sym_mul(A, v) = sympy_mat(A) * sympy_vec(v)
sym_add(A, B) = sympy_mat(A) + sympy_mat(B)
sym_pow(A, k) = sympy_mat(A) ^ k

sym_is_zero(x) = PythonCall.pyconvert(Bool, _sympy_module().simplify(x).is_zero_matrix)
sym_eq(A, B) = sym_is_zero(sympy_mat(A) - sympy_mat(B))

sym_vec_zero(v) = all(PythonCall.pyconvert(Bool, _sympy_module().simplify(e) == 0) for e in v)

sympy_to_julia_vec(x) = x isa PythonCall.Py ? Base.invokelatest(PythonCall.pyconvert, Vector{Any}, x) : x
sympy_to_julia_mat(x) = x isa PythonCall.Py ? Base.invokelatest(PythonCall.pyconvert, Matrix{Any}, x) : x

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
    sympy_subs_numeric(A, subs) -> Union{Py, AbstractArray}

Substitute `subs` into a SymPy matrix or Julia matrix convertible to SymPy.
`subs` can be a Dict or list of pairs. Returns a SymPy matrix if symbols remain;
otherwise returns a Julia numeric array.
"""
function sympy_subs_numeric(A, subs)
    sympy = _sympy_module()
    symA = sympy_mat(A)
    sub_list = subs isa AbstractDict ? collect(pairs(subs)) : subs
    subbed = symA.subs(sub_list)
    free = Base.invokelatest(PythonCall.pygetattr, subbed, "free_symbols")
    blen = Base.invokelatest(PythonCall.pybuiltins.len, free)
    nfree = Base.invokelatest(PythonCall.pyconvert, Int, blen)
    if nfree != 0
        return subbed
    end
    M = sympy_to_julia_mat(subbed)
    try
        num = map(_sympy_scalar_to_julia, M)
        return _promote_matrix(num)
    catch
        return subbed
    end
end

end
