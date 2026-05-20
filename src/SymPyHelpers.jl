module SymPyHelpers

export sym_mat, sym_vec, sym_zero, sym_mul, sym_add, sym_pow, sym_eq, sym_is_zero, sym_vec_zero
export sym_to_julia_vec, sym_to_julia_mat, sym_subs_numeric

using ..GenLAProblems: ensure_pythoncall!

const _sympy = Ref{Any}(nothing)

_py() = ensure_pythoncall!()
_is_py(x) = x isa getfield(_py(), :Py)

function _sympy_module()
    if _sympy[] === nothing
        _sympy[] = getfield(parentmodule(@__MODULE__), :sympy)
    end
    return _sympy[]
end

function _to_sympy_entry(sympy, x)
    if _is_py(x)
        return x
    elseif x isa Rational
        return sympy.Rational(numerator(x), denominator(x))
    elseif x isa Complex{<:Rational}
        r = sympy.Rational(numerator(real(x)), denominator(real(x)))
        i = sympy.Rational(numerator(imag(x)), denominator(imag(x)))
        return r + i * sympy.I
    elseif x isa Complex
        return real(x) + imag(x) * sympy.I
    else
        return x
    end
end

function _sympy_matrix_from_array(x::AbstractArray)
    py = _py()
    sympy = _sympy_module()
    # Normalize vectors to a column matrix to avoid Matrix(::Vector{Py}) errors.
    mat = x isa AbstractVector ? reshape(x, :, 1) : Matrix(x)
    rows = Vector{Any}(undef, size(mat, 1))
    for i in 1:size(mat, 1)
        row = Vector{Any}(undef, size(mat, 2))
        for j in 1:size(mat, 2)
            row[j] = _to_sympy_entry(sympy, mat[i, j])
        end
        rows[i] = row
    end
    # Build explicit Python lists to avoid NULL conversions in nested Julia arrays.
    pyrows = Base.invokelatest(py.pylist, [
        Base.invokelatest(py.pylist, r) for r in rows
    ])
    return sympy.Matrix(pyrows)
end

sym_mat(x) = _is_py(x) ? x : (x isa AbstractArray ? _sympy_matrix_from_array(x) : _sympy_module().Matrix(x))
sym_vec(x) = _is_py(x) ? x : (x isa AbstractArray ? _sympy_matrix_from_array(x) : _sympy_module().Matrix(x))
sym_zero() = _sympy_module().Integer(0)

sym_mul(A, v) = sym_mat(A) * sym_vec(v)
sym_add(A, B) = sym_mat(A) + sym_mat(B)
sym_pow(A, k) = sym_mat(A) ^ k

function sym_is_zero(x)
    py = _py()
    return Base.invokelatest(py.pyconvert, Bool, _sympy_module().simplify(x).is_zero_matrix)
end
sym_eq(A, B) = sym_is_zero(sym_mat(A) - sym_mat(B))

function sym_vec_zero(v)
    py = _py()
    return all(Base.invokelatest(py.pyconvert, Bool, _sympy_module().simplify(e) == 0) for e in v)
end

function sym_to_julia_vec(x)
    py = _py()
    return _is_py(x) ? Base.invokelatest(py.pyconvert, Vector{Any}, x) : x
end

function sym_to_julia_mat(x)
    py = _py()
    return _is_py(x) ? Base.invokelatest(py.pyconvert, Matrix{Any}, x) : x
end

function _sympy_scalar_to_julia(x)
    if !_is_py(x)
        return x
    end
    py = _py()
    sympy = _sympy_module()
    is_int = Base.invokelatest(py.pyconvert, Bool, Base.invokelatest(py.pygetattr, x, "is_Integer"))
    if is_int
        return Base.invokelatest(py.pyconvert, Int, x)
    end
    is_rat = Base.invokelatest(py.pyconvert, Bool, Base.invokelatest(py.pygetattr, x, "is_Rational"))
    if is_rat
        p = Base.invokelatest(py.pyconvert, Int, Base.invokelatest(py.pygetattr, x, "p"))
        q = Base.invokelatest(py.pyconvert, Int, Base.invokelatest(py.pygetattr, x, "q"))
        return p//q
    end
    return Base.invokelatest(py.pyconvert, Float64, sympy.N(x))
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
    py = _py()
    symA = sym_mat(A)
    sub_list = subs isa AbstractDict ? collect(pairs(subs)) : subs
    if sub_list isa Pair
        sub_list = [sub_list]
    end
    if sub_list isa AbstractVector
        sub_list = [(p.first, p.second) for p in sub_list]
    end
    subbed = symA.subs(sub_list)
    free = Base.invokelatest(py.pygetattr, subbed, "free_symbols")
    blen = Base.invokelatest(py.pybuiltins.len, free)
    nfree = Base.invokelatest(py.pyconvert, Int, blen)
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
