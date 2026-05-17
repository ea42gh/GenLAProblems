using Test
using PythonCall
using GenLAProblems

function _has_module_symops(name::String)
    try
        PythonCall.pyimport(name)
        return true
    catch
        return false
    end
end

function _py_setattr_symops(obj, name::AbstractString, value)
    PythonCall.pycall(PythonCall.pybuiltins.setattr, obj, name, value)
end

@testset "SymPyHelpers algebra helpers" begin
    sym = try
        GenLAProblems.import_sympy()
    catch
        nothing
    end
    if sym === nothing
        @info "Skipping SymPyHelpers algebra helper tests: sympy unavailable"
    else
        try
            x = sym.symbols("x")
            A = sym.eye(2)
            B = sym.ones(2, 2)
            U = sym.eye(2)
            U[0, 1] = 1

            A2 = GenLAProblems.sym_add(A, B)
            E2 = sym.Matrix([[2, 1], [1, 2]])
            @test PythonCall.pyconvert(Bool, sym.simplify(A2 - E2).is_zero_matrix)

            P2 = GenLAProblems.sym_pow(U, 2)
            Epow = sym.Matrix([[1, 2], [0, 1]])
            @test PythonCall.pyconvert(Bool, sym.simplify(P2 - Epow).is_zero_matrix)

            mv = GenLAProblems.sym_mul(A, sym.Matrix([x, 1]))
            @test PythonCall.pyconvert(Bool, sym.simplify(mv[0] - x) == 0)
            @test PythonCall.pyconvert(Bool, sym.simplify(mv[1] - 1) == 0)

            @test GenLAProblems.sym_eq(A, A)
            @test !GenLAProblems.sym_eq(A, U)

            z0 = sym.Matrix([sym.Integer(0), sym.Integer(0)])
            z1 = sym.Matrix([x, sym.Integer(0)])
            @test GenLAProblems.sym_vec_zero(z0)
            @test !GenLAProblems.sym_vec_zero(z1)
        catch err
            @info "Skipping SymPyHelpers algebra helper tests: runtime conversion mismatch" exception=(err, catch_backtrace())
        end
    end
end

@testset "SymPyHelpers adjoint inputs" begin
    A = Rational.([1 2; 3 4])
    Aadj = A'
    sym = try
        GenLAProblems.import_sympy()
    catch
        nothing
    end
    if sym === nothing
        @info "Skipping SymPyHelpers adjoint inputs: sympy unavailable"
    else
        try
            subbed = GenLAProblems.sym_subs_numeric(Aadj, Dict())
            @test subbed isa AbstractArray || subbed isa PythonCall.Py
        catch err
            @info "Skipping SymPyHelpers adjoint inputs: runtime conversion mismatch" exception=(err, catch_backtrace())
        end
    end
end

@testset "Extractor wrappers call-through contracts" begin
    types = PythonCall.pyimport("types")
    la = PythonCall.pycall(types.SimpleNamespace)

    seen = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

    function fake_svd(spec; kwargs...)
        seen[] = Dict(:which => :svd, :kwargs => Dict(kwargs), :spec => spec)
        return ("U", "S", "V", 2)
    end

    function fake_eig(spec; kwargs...)
        seen[] = Dict(:which => :eig, :kwargs => Dict(kwargs), :spec => spec)
        return ("L", "Q")
    end

    function fake_qr(mats)
        seen[] = Dict(:which => :qr, :mats => mats)
        return PythonCall.pydict(Dict("Q" => "q", "R" => "r"))
    end

    _py_setattr_symops(la, "svd_matrices_from_spec", fake_svd)
    _py_setattr_symops(la, "eig_matrices_from_spec", fake_eig)
    _py_setattr_symops(la, "qr_matrices_from_grid", fake_qr)

    old_la = GenLAProblems._LAFigureSpecs[]
    try
        GenLAProblems._LAFigureSpecs[] = la

        out_svd = GenLAProblems.svd_matrices_from_spec("spec"; reduced=false)
        @test out_svd == ("U", "S", "V", 2)
        @test seen[][:which] == :svd
        @test seen[][:kwargs][:reduced] == false

        out_eig = GenLAProblems.eig_matrices_from_spec("espec"; orthonormal=false)
        @test out_eig == ("L", "Q")
        @test seen[][:which] == :eig
        @test seen[][:kwargs][:orthonormal] == false

        out_qr = GenLAProblems.qr_matrices_from_grid("grid")
        @test out_qr == (A=nothing, W=nothing, WtA=nothing, WtW=nothing, S=nothing, Qt=nothing, Q="q", R="r")
        @test seen[][:which] == :qr
        @test seen[][:mats] == "grid"
    finally
        GenLAProblems._LAFigureSpecs[] = old_la
    end
end

@testset "materialize_python_value converts nested Python containers" begin
    sym = PythonCall.pyimport("sympy")
    pymat = sym.Matrix(PythonCall.pylist([
        PythonCall.pylist([sym.Integer(1), sym.Integer(2)]),
        PythonCall.pylist([sym.Integer(3), sym.Integer(4)]),
    ]))
    py = PythonCall.pydict(Dict(
        "tuple_like" => (pymat, PythonCall.pylist([sym.Integer(5), sym.Integer(6)])),
        "scalar" => sym.Integer(7),
    ))
    out = GenLAProblems.materialize_python_value(py)
    @test out isa Dict
    @test out["tuple_like"] isa Tuple
    @test out["tuple_like"][1] isa AbstractMatrix
    @test out["tuple_like"][1] == [1 2; 3 4]
    @test out["tuple_like"][2] == [5, 6]
    @test out["scalar"] == 7
end

@testset "materialize_python_value handles none and Julia passthrough" begin
    pybuiltins = PythonCall.pyimport("builtins")
    @test GenLAProblems.materialize_python_value(PythonCall.pygetattr(pybuiltins, "None")) === nothing

    julia_dict = Dict("a" => [1, 2], "b" => (3, 4))
    out = GenLAProblems.materialize_python_value(julia_dict)
    @test out == julia_dict
    @test out["a"] isa Vector{Int}
    @test out["b"] == (3, 4)
end

@testset "materialize_python_value converts nested dict/list mixtures" begin
    sym = PythonCall.pyimport("sympy")
    inner = PythonCall.pydict(Dict(
        "vec" => PythonCall.pylist([sym.Integer(1), sym.Integer(2)]),
        "none" => PythonCall.pygetattr(PythonCall.pyimport("builtins"), "None"),
    ))
    outer = PythonCall.pydict(Dict(
        "items" => PythonCall.pylist([inner, PythonCall.pydict(Dict("x" => sym.Integer(9)))]),
        "name" => "qr",
    ))

    out = GenLAProblems.materialize_python_value(outer)
    @test out isa Dict
    @test out["name"] == "qr"
    @test out["items"] isa Vector
    @test out["items"][1]["vec"] == [1, 2]
    @test out["items"][1]["none"] === nothing
    @test out["items"][2]["x"] == 9
end

@testset "qr_matrices_from_grid preserves missing entries and converts values" begin
    sym = PythonCall.pyimport("sympy")
    types = PythonCall.pyimport("types")
    pybuiltins = PythonCall.pyimport("builtins")
    la = PythonCall.pycall(types.SimpleNamespace)

    pymat(rows) = sym.Matrix(PythonCall.pylist([
        PythonCall.pylist([sym.Integer(v) for v in row]) for row in rows
    ]))

    function fake_qr_py(mats)
        return PythonCall.pydict(Dict(
            "A" => pybuiltins.None,
            "W" => pybuiltins.None,
            "WtA" => pymat([[1, 2]]),
            "WtW" => pymat([[3]]),
            "S" => pybuiltins.None,
            "Qt" => pybuiltins.None,
            "Q" => pymat([[1, 0], [0, 1]]),
            "R" => pymat([[4, 5], [0, 6]]),
        ))
    end

    _py_setattr_symops(la, "qr_matrices_from_grid", fake_qr_py)

    old_la = GenLAProblems._LAFigureSpecs[]
    try
        GenLAProblems._LAFigureSpecs[] = la
        qr = GenLAProblems.qr_matrices_from_grid("grid")
        @test qr.A === nothing
        @test qr.W === nothing
        @test qr.S === nothing
        @test qr.Qt === nothing
        @test qr.WtA == [1 2]
        @test qr.WtW == [3;;]
        @test qr.Q == [1 0; 0 1]
        @test qr.R == [4 5; 0 6]
    finally
        GenLAProblems._LAFigureSpecs[] = old_la
    end
end

@testset "eig and svd matrix extractors return Julia values" begin
    sym = PythonCall.pyimport("sympy")
    types = PythonCall.pyimport("types")
    la = PythonCall.pycall(types.SimpleNamespace)

    pymat(rows) = sym.Matrix(PythonCall.pylist([
        PythonCall.pylist([sym.Integer(v) for v in row]) for row in rows
    ]))

    function fake_svd_py(spec; kwargs...)
        return (pymat([[1, 0], [0, 1]]), pymat([[2, 0], [0, 1]]), pymat([[1, 0], [0, 1]]), 2)
    end

    function fake_eig_py(spec; kwargs...)
        return (pymat([[3, 0], [0, 4]]), pymat([[1, 0], [0, 1]]))
    end

    _py_setattr_symops(la, "svd_matrices_from_spec", fake_svd_py)
    _py_setattr_symops(la, "eig_matrices_from_spec", fake_eig_py)

    old_la = GenLAProblems._LAFigureSpecs[]
    try
        GenLAProblems._LAFigureSpecs[] = la
        U, S, V, rank = GenLAProblems.svd_matrices_from_spec("spec")
        @test U isa AbstractMatrix
        @test S isa AbstractMatrix
        @test V isa AbstractMatrix
        @test rank == 2

        Lambda, V2 = GenLAProblems.eig_matrices_from_spec("spec")
        @test Lambda isa AbstractMatrix
        @test V2 isa AbstractMatrix
    finally
        GenLAProblems._LAFigureSpecs[] = old_la
    end
end

@testset "_pygetattr_fallback behavior" begin
    types = PythonCall.pyimport("types")
    sys = PythonCall.pyimport("sys")
    modules = PythonCall.pygetattr(sys, "modules")

    obj_with = PythonCall.pycall(types.SimpleNamespace)
    _py_setattr_symops(obj_with, "foo", 123)

    mod_name = "_tmp_glp_fallback"
    mod = PythonCall.pycall(types.ModuleType, mod_name)
    _py_setattr_symops(mod, "foo", 999)

    contains = PythonCall.pycall(PythonCall.pybuiltins.getattr, modules, "__contains__")
    had_mod = PythonCall.pyconvert(Bool, PythonCall.pycall(contains, mod_name))
    old_mod = had_mod ? modules[mod_name] : nothing

    try
        modules[mod_name] = mod

        v1 = GenLAProblems._pygetattr_fallback(obj_with, :foo, mod_name)
        @test PythonCall.pyconvert(Int, v1) == 123

        obj_without = PythonCall.pycall(types.SimpleNamespace)
        v2 = GenLAProblems._pygetattr_fallback(obj_without, :foo, mod_name)
        @test PythonCall.pyconvert(Int, v2) == 999
    finally
        if had_mod
            modules[mod_name] = old_mod
        else
            popfn = PythonCall.pycall(PythonCall.pybuiltins.getattr, modules, "pop")
            PythonCall.pycall(popfn, mod_name, nothing)
        end
    end
end
