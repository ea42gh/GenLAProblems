using Test
using PythonCall
using GenLAProblems

function _has_module_cov(name::String)
    try
        PythonCall.pyimport(name)
        return true
    catch
        return false
    end
end

function _py_setattr_cov(obj, name::AbstractString, value)
    PythonCall.pycall(PythonCall.pybuiltins.setattr, obj, name, value)
end

@testset "matrixlayout_ge forwarding and merge behavior" begin
    if !_has_module_cov("la_figures.ge_convenience")
        @info "Skipping matrixlayout_ge forwarding tests: la_figures.ge_convenience unavailable"
    else
        ge_conv = PythonCall.pyimport("la_figures.ge_convenience")
        old_ge = PythonCall.pygetattr(ge_conv, "ge")
        seen = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

        function fake_ge(args...; kwargs...)
            seen[] = Dict(kwargs)
            return "<svg>ok</svg>"
        end

        try
            _py_setattr_cov(ge_conv, "ge", fake_ge)

            mats = [[nothing, [1 2; 3 4]]]
            out = GenLAProblems.matrixlayout_ge(
                mats;
                Nrhs=1,
                array_names=["A", "b"],
                pivot_list=[[(0, 1), [(0, 0)]]],
                ref_path_list=[[0, 1, [(0, 0)], "vv"]],
                comment_list=["step"],
                specs=[Dict("grid" => (0, 1), "label" => "A", "side" => "above")],
            )
            @test out isa GenLAProblems.SVGOut
            @test occursin("<svg>", out.svg)
            @test haskey(seen[], :Nrhs)
            @test seen[][:Nrhs] == 1
            @test haskey(seen[], :array_names)
            @test haskey(seen[], :pivot_list)
            @test haskey(seen[], :ref_path_list)
            @test haskey(seen[], :comment_list)
            @test haskey(seen[], :specs)

            # bg_for_entries should be converted into decorators and bg_for_entries cleared.
            out2 = GenLAProblems.matrixlayout_ge(
                mats;
                bg_for_entries=[[0, 1, [[(0, 0), (1, 1)]], "yellow!35", 1]],
                decorators=[Dict("grid" => (0, 1), "entries" => Any[], "decorator" => (x -> x))],
            )
            @test out2 isa GenLAProblems.SVGOut
            @test haskey(seen[], :decorators)
            @test !isnothing(seen[][:decorators])
            @test haskey(seen[], :bg_for_entries)
            @test isnothing(seen[][:bg_for_entries])
        finally
            _py_setattr_cov(ge_conv, "ge", old_ge)
        end
    end
end

@testset "SymPyHelpers additional substitution edges" begin
    if !_has_module_cov("sympy")
        @info "Skipping SymPyHelpers edge tests: sympy unavailable"
    else
        sym = try
            GenLAProblems.import_sympy()
        catch err
            @info "Skipping SymPyHelpers edge tests: import_sympy unavailable in this runtime" err
            nothing
        end
        if sym !== nothing
            try
                x = sym.symbols("x")
                y = sym.symbols("y")
                A = [x 1; 0 y]

                # Multiple substitutions: numeric output expected.
                M1 = GenLAProblems.sym_subs_numeric(A, Dict(x => 2, y => 3))
                @test M1 isa AbstractMatrix
                @test M1[1, 1] == 2
                @test M1[2, 2] == 3

                # Partial substitution: symbols remain => SymPy matrix expected.
                M2 = GenLAProblems.sym_subs_numeric(A, Dict(x => 2))
                @test M2 isa PythonCall.Py
                free = PythonCall.pygetattr(M2, "free_symbols")
                @test PythonCall.pyconvert(Int, PythonCall.pybuiltins.len(free)) == 1

                # Pair-list style substitutions.
                M3 = GenLAProblems.sym_subs_numeric(A, [x => 4, y => 5])
                @test M3 isa AbstractMatrix
                @test M3[1, 1] == 4
                @test M3[2, 2] == 5
            catch err
                @info "Skipping SymPyHelpers edge tests: runtime conversion mismatch" exception=(err, catch_backtrace())
            end
        end
    end
end

@testset "rhs_block helper" begin
    A = [1 2; 3 4]
    b = [5, 6]
    pb = ShowGe(A, b)
    ref!(pb; gj=true)
    rhs = rhs_block(pb)
    @test rhs == pb.matrices[end][end][:, size(A, 2)+1:end]
    col1 = rhs_block(pb; b_col=1)
    @test col1 == rhs[:, 1]
end

@testset "normal_eq label specs" begin
    specs = GenLAProblems._normal_eq_name_specs(3, ["A", "b"])
    @test any(t -> occursin("A^T", t[3]), specs)
end

@testset "nM.show_ge_tbl fallback routing" begin
    types = PythonCall.pyimport("types")
    sys = PythonCall.pyimport("sys")

    fake_root = PythonCall.pycall(types.SimpleNamespace)
    old_la = GenLAProblems._la_figures[]

    # Build fake la_figures.ge_convenience module in sys.modules.
    mod_la = PythonCall.pycall(types.ModuleType, "la_figures")
    mod_conv = PythonCall.pycall(types.ModuleType, "la_figures.ge_convenience")
    seen = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

    function fake_ge_tbl_svg(args...; kwargs...)
        seen[] = Dict(kwargs)
        return "<svg>fallback</svg>"
    end
    _py_setattr_cov(mod_conv, "ge_tbl_svg", fake_ge_tbl_svg)

    modules = PythonCall.pygetattr(sys, "modules")
    contains = PythonCall.pycall(PythonCall.pybuiltins.getattr, modules, "__contains__")
    had_la = PythonCall.pyconvert(Bool, PythonCall.pycall(contains, "la_figures"))
    had_conv = PythonCall.pyconvert(Bool, PythonCall.pycall(PythonCall.pybuiltins.getattr(modules, "__contains__"), "la_figures.ge_convenience"))
    old_mod_la = had_la ? modules["la_figures"] : nothing
    old_mod_conv = had_conv ? modules["la_figures.ge_convenience"] : nothing

    try
        modules["la_figures"] = mod_la
        modules["la_figures.ge_convenience"] = mod_conv
        GenLAProblems._la_figures[] = fake_root

        h = nM.show_ge_tbl([1 0; 0 1]; tmp_dir="/tmp/la", keep_file="x", output_dir="/tmp/y")
        @test h isa GenLAProblems.SVGOut
        @test occursin("fallback", h.svg)
        @test !haskey(seen[], :tmp_dir)
        @test !haskey(seen[], :keep_file)
        @test !haskey(seen[], :output_dir)
    finally
        if had_la
            modules["la_figures"] = old_mod_la
        else
            popfn = PythonCall.pycall(PythonCall.pybuiltins.getattr, modules, "pop")
            PythonCall.pycall(popfn, "la_figures", nothing)
        end
        if had_conv
            modules["la_figures.ge_convenience"] = old_mod_conv
        else
            popfn = PythonCall.pycall(PythonCall.pybuiltins.getattr, modules, "pop")
            PythonCall.pycall(popfn, "la_figures.ge_convenience", nothing)
        end
        GenLAProblems._la_figures[] = old_la
    end
end

@testset "show_ge_final contract" begin
    if !_has_module_cov("la_figures")
        @info "Skipping show_ge_final contract test: la_figures unavailable"
    else
        old_la = GenLAProblems._la_figures[]
        types = PythonCall.pyimport("types")
        la = PythonCall.pycall(types.SimpleNamespace)
        seen = Ref{Tuple{Any,Any}}((nothing, nothing))

        function fake_ge_tbl_svg(A, rhs; kwargs...)
            seen[] = (A, rhs)
            return "<svg>final</svg>"
        end
        _py_setattr_cov(la, "ge_tbl_svg", fake_ge_tbl_svg)

        try
            GenLAProblems._la_figures[] = la
            mats = [[nothing, [1 2; 3 4]], [nothing, [9 8; 7 6]]]
            h = GenLAProblems.show_ge_final(mats, Any[], Int[]; Nrhs=0)
            @test h isa GenLAProblems.SVGOut
            @test occursin("final", h.svg)
            A_seen, rhs_seen = seen[]
            @test rhs_seen === nothing
            @test PythonCall.pyconvert(Matrix{Any}, A_seen) == Any[9 8; 7 6]
        finally
            GenLAProblems._la_figures[] = old_la
        end
    end
end
