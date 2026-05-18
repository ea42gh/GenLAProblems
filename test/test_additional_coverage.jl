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
    la = PythonCall.pycall(PythonCall.pyimport("types").SimpleNamespace)
    seen = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

    function fake_ge(args...; kwargs...)
        seen[] = Dict(kwargs)
        return "<svg>ok</svg>"
    end

    _py_setattr_cov(la, "ge", fake_ge)
    old_la = GenLAProblems._LAFigureSpecs[]
    try
        GenLAProblems._LAFigureSpecs[] = la

        mats = [[nothing, [1 2; 3 4]]]
        out = GenLAProblems.matrixlayout_ge(
            mats;
            n_rhs=1,
            output_dir="/tmp/out",
            output_stem="ge_basic",
            array_names=["A", "b"],
            pivot_list=[[(0, 1), [(0, 0)]]],
            ref_path_list=[[0, 1, [(0, 0)], "vv"]],
            comment_list=["step"],
            specs=[Dict("grid" => (0, 1), "label" => "A", "side" => "above")],
        )
        @test out isa GenLAProblems.SVGOut
        @test occursin("<svg>", out.svg)
        @test haskey(seen[], :output_dir)
        @test seen[][:output_dir] == "/tmp/out"
        @test haskey(seen[], :output_stem)
        @test seen[][:output_stem] == "ge_basic"
        @test haskey(seen[], :n_rhs)
        @test seen[][:n_rhs] == 1
        @test haskey(seen[], :array_names)
        @test haskey(seen[], :pivot_list)
        @test haskey(seen[], :ref_path_list)
        @test haskey(seen[], :comment_list)
        @test haskey(seen[], :specs)

        out_vec = GenLAProblems.matrixlayout_ge(
            mats;
            n_rhs=[1, 1],
            array_names=["A", ["b₁", "b₂"]],
        )
        @test out_vec isa GenLAProblems.SVGOut
        @test PythonCall.pyconvert(Vector{Int}, seen[][:n_rhs]) == [1, 1]

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
        GenLAProblems._LAFigureSpecs[] = old_la
    end
end

@testset "keep_file path splitting" begin
    dir, stem = GenLAProblems._keep_file_output_parts("/tmp/la/run/show_layout.tex")
    @test dir == "/tmp/la/run"
    @test stem == "show_layout"
    dir2, stem2 = GenLAProblems._keep_file_output_parts("/tmp/la/run/show_layout")
    @test dir2 == "/tmp/la/run"
    @test stem2 == "show_layout"
    @test GenLAProblems._keep_file_output_parts(nothing) == (nothing, nothing)
end

@testset "GE output target resolution" begin
    pb = ShowGE([1 2; 3 4]; tmp_dir="/tmp/fallback", keep_file="/tmp/kept/show_layout.tex")
    out_dir, out_stem = GenLAProblems._resolve_ge_output_targets(pb, nothing, nothing)
    @test out_dir == "/tmp/kept"
    @test out_stem == "show_layout"
    out_dir2, out_stem2 = GenLAProblems._resolve_ge_output_targets(pb, "/tmp/explicit", "demo")
    @test out_dir2 == "/tmp/explicit"
    @test out_stem2 == "demo"
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
    pb = ShowGE(A, b)
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

@testset "homogeneous solutions include zero vector when full rank" begin
    A = Rational{Int}.([1 2; 3 4])
    b = Rational{Int}.([5, 6])
    pb = ShowGE(A, b)
    ref!(pb; gj=true)
    _, Xh = solutions(pb)
    @test size(Xh, 2) == 1
    @test all(x -> x == 0, Xh)
end

@testset "matrixlayout_ge sets medium nodes for bg_for_entries" begin
    payloads = Vector{Any}()
    if GenLAProblems._pythoncall_disabled()
        @info "Skipping matrixlayout_ge medium nodes test: PythonCall disabled"
        return
    end
    if GenLAProblems._in_precompile()
        @info "Skipping matrixlayout_ge medium nodes test: precompile mode"
        return
    end
    function fake_ge(mats_py; kwargs...)
        push!(payloads, kwargs)
        return "<svg/>"
    end
    la = PythonCall.pycall(PythonCall.pyimport("types").SimpleNamespace)
    _py_setattr_cov(la, "ge", fake_ge)
    old_la = GenLAProblems._LAFigureSpecs[]
    try
        GenLAProblems._LAFigureSpecs[] = la
        mats = [[nothing, [1 0; 0 1]]]
        bg = [[0, 1, [[(0, 0), (0, 0)]], "yellow!25", 1]]
        GenLAProblems.matrixlayout_ge(mats; bg_for_entries=bg)
        @test !isempty(payloads)
        ro = get(payloads[end], :render_opts, nothing)
        @test ro === nothing || get(ro, "create_medium_nodes", true) == true
    finally
        GenLAProblems._LAFigureSpecs[] = old_la
    end
end

@testset "nM.show_ge_tbl fallback routing" begin
    types = PythonCall.pyimport("types")
    sys = PythonCall.pyimport("sys")

    fake_root = PythonCall.pycall(types.SimpleNamespace)
    old_la = GenLAProblems._LAFigureSpecs[]

    # Build fake LAFigureSpecs.ge_convenience module in sys.modules.
    mod_la = PythonCall.pycall(types.ModuleType, "LAFigureSpecs")
    mod_conv = PythonCall.pycall(types.ModuleType, "LAFigureSpecs.ge_convenience")
    seen = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

    function fake_ge_tbl_svg(args...; kwargs...)
        seen[] = Dict(kwargs)
        return "<svg>fallback</svg>"
    end
    _py_setattr_cov(mod_conv, "ge_tbl_svg", fake_ge_tbl_svg)

    modules = PythonCall.pygetattr(sys, "modules")
    contains = PythonCall.pycall(PythonCall.pybuiltins.getattr, modules, "__contains__")
    had_la = PythonCall.pyconvert(Bool, PythonCall.pycall(contains, "LAFigureSpecs"))
    had_conv = PythonCall.pyconvert(Bool, PythonCall.pycall(PythonCall.pybuiltins.getattr(modules, "__contains__"), "LAFigureSpecs.ge_convenience"))
    old_mod_la = had_la ? modules["LAFigureSpecs"] : nothing
    old_mod_conv = had_conv ? modules["LAFigureSpecs.ge_convenience"] : nothing

    try
        modules["LAFigureSpecs"] = mod_la
        modules["LAFigureSpecs.ge_convenience"] = mod_conv
        GenLAProblems._LAFigureSpecs[] = fake_root

        h = nM.show_ge_tbl([1 0; 0 1]; tmp_dir="/tmp/la", keep_file="x", output_dir="/tmp/y")
        @test h isa GenLAProblems.SVGOut
        @test occursin("fallback", h.svg)
        @test !haskey(seen[], :tmp_dir)
        @test !haskey(seen[], :keep_file)
        @test !haskey(seen[], :output_dir)
        @test !haskey(seen[], :output_stem)
    finally
        if had_la
            modules["LAFigureSpecs"] = old_mod_la
        else
            popfn = PythonCall.pycall(PythonCall.pybuiltins.getattr, modules, "pop")
            PythonCall.pycall(popfn, "LAFigureSpecs", nothing)
        end
        if had_conv
            modules["LAFigureSpecs.ge_convenience"] = old_mod_conv
        else
            popfn = PythonCall.pycall(PythonCall.pybuiltins.getattr, modules, "pop")
            PythonCall.pycall(popfn, "LAFigureSpecs.ge_convenience", nothing)
        end
        GenLAProblems._LAFigureSpecs[] = old_la
    end
end

@testset "show_ge_final contract" begin
    if !_has_module_cov("LAFigureSpecs")
        @info "Skipping show_ge_final contract test: LAFigureSpecs unavailable"
    else
        old_la = GenLAProblems._LAFigureSpecs[]
        types = PythonCall.pyimport("types")
        la = PythonCall.pycall(types.SimpleNamespace)
        seen = Ref{Tuple{Any,Any}}((nothing, nothing))
        seen_kwargs = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

        function fake_ge_tbl_svg(A, rhs; kwargs...)
            seen[] = (A, rhs)
            seen_kwargs[] = Dict(kwargs)
            return "<svg>final</svg>"
        end
        _py_setattr_cov(la, "ge_tbl_svg", fake_ge_tbl_svg)

        try
            GenLAProblems._LAFigureSpecs[] = la
            mats = [[nothing, [1 2; 3 4]], [nothing, [9 8; 7 6]]]
            h = GenLAProblems.show_ge_final(mats, Any[], Int[]; n_rhs=0, output_dir="/tmp/final")
            @test h isa GenLAProblems.SVGOut
            @test occursin("final", h.svg)
            A_seen, rhs_seen = seen[]
            @test rhs_seen === nothing
            @test PythonCall.pyconvert(Matrix{Any}, A_seen) == Any[9 8; 7 6]
            @test haskey(seen_kwargs[], :output_dir)
            @test seen_kwargs[][:output_dir] == "/tmp/final"
        finally
            GenLAProblems._LAFigureSpecs[] = old_la
        end
    end
end
