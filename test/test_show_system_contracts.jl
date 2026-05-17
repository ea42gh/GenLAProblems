using Test
using PythonCall
using GenLAProblems

function _py_setattr_sys(obj, name::AbstractString, value)
    PythonCall.pycall(PythonCall.pybuiltins.setattr, obj, name, value)
end

@testset "show_system contracts (mocked la_figures + matrixlayout.backsubst)" begin
    types = PythonCall.pyimport("types")
    sys = PythonCall.pyimport("sys")
    modules = PythonCall.pygetattr(sys, "modules")

    # Fake la_figures root (only linear_system_tex needed).
    la = PythonCall.pycall(types.SimpleNamespace)
    seen_tex = Ref{NamedTuple{(:A,:b,:var_name),Tuple{Any,Any,Any}}}((A=nothing, b=nothing, var_name=nothing))
    seen_svg = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

    function fake_linear_system_tex(A, b; var_name="x")
        seen_tex[] = (A=A, b=b, var_name=var_name)
        return "\\systeme{x_1=1}"
    end
    _py_setattr_sys(la, "linear_system_tex", fake_linear_system_tex)

    # Fake matrixlayout.backsubst module.
    mod_bs = PythonCall.pycall(types.ModuleType, "matrixlayout.backsubst")
    function fake_backsubst_svg(; kwargs...)
        seen_svg[] = Dict(kwargs)
        return "<svg>system</svg>"
    end
    _py_setattr_sys(mod_bs, "backsubst_svg", fake_backsubst_svg)

    contains = PythonCall.pycall(PythonCall.pybuiltins.getattr, modules, "__contains__")
    had_ml = PythonCall.pyconvert(Bool, PythonCall.pycall(contains, "matrixlayout"))
    old_ml = had_ml ? modules["matrixlayout"] : nothing
    had_bs = PythonCall.pyconvert(Bool, PythonCall.pycall(contains, "matrixlayout.backsubst"))
    old_bs = had_bs ? modules["matrixlayout.backsubst"] : nothing
    old_la = GenLAProblems._la_figures[]

    try
        modules["matrixlayout"] = PythonCall.pycall(types.ModuleType, "matrixlayout")
        modules["matrixlayout.backsubst"] = mod_bs
        GenLAProblems._la_figures[] = la

        # Integer path: use selected RHS column.
        A = [1 2; 3 4]
        B = [10 11; 20 21]
        pb = ShowGE(A, B; output_dir="/tmp/la")
        h = show_system(pb; b_col=2, var_name="z", fig_scale=1.7, output_dir="/tmp/sys")
        @test h isa GenLAProblems.SVGOut
        @test occursin("system", h.svg)
        @test seen_tex[].var_name == "z"
        @test PythonCall.pyconvert(Matrix{Any}, seen_tex[].A) == Any[1 2; 3 4]
        @test PythonCall.pyconvert(Vector{Any}, seen_tex[].b) == Any[11, 21]
        @test seen_svg[][:show_system] == true
        @test seen_svg[][:show_cascade] == false
        @test seen_svg[][:show_solution] == false
        @test seen_svg[][:fig_scale] == 1.7
        @test seen_svg[][:output_dir] == "/tmp/sys"

        # Rational path: A and b are encoded as tuple entries.
        Ar = Rational{Int}.(A)
        Br = Rational{Int}.(B)
        pbr = ShowGE{Rational{Int}}(A, B; output_dir="/tmp/la")
        _ = show_system(pbr; b_col=1, var_name="x", output_dir="/tmp/rational")
        A_r = PythonCall.pyconvert(Matrix{Any}, seen_tex[].A)
        b_r = PythonCall.pyconvert(Vector{Any}, seen_tex[].b)
        @test A_r[1, 1] == (1, 1)
        @test b_r[1] == (10, 1)

        # Complex{Rational} path: should still pass tuple-like encoded entries and render.
        Ac = Complex{Int}.(A)
        Bc = Complex{Int}.(B)
        pbc = ShowGE{Complex{Rational{Int}}}(Ac, Bc; output_dir="/tmp/la")
        h3 = show_system(pbc; b_col=1, var_name="w", output_dir="/tmp/complex")
        @test h3 isa GenLAProblems.SVGOut
        @test seen_tex[].var_name == "w"
    finally
        if had_bs
            modules["matrixlayout.backsubst"] = old_bs
        else
            popfn = PythonCall.pycall(PythonCall.pybuiltins.getattr, modules, "pop")
            PythonCall.pycall(popfn, "matrixlayout.backsubst", nothing)
        end
        if had_ml
            modules["matrixlayout"] = old_ml
        else
            popfn = PythonCall.pycall(PythonCall.pybuiltins.getattr, modules, "pop")
            PythonCall.pycall(popfn, "matrixlayout", nothing)
        end
        GenLAProblems._la_figures[] = old_la
    end
end
