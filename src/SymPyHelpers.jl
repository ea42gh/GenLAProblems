module SymPyHelpers

export sympy_mat, sympy_vec, sympy_zero, sym_mul, sym_add, sym_pow, sym_eq, sym_is_zero, sym_vec_zero
export sympy_to_julia_vec, sympy_to_julia_mat

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

end
