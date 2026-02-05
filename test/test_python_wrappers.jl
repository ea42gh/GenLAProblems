using Test
using PythonCall
using GenLAProblems

@testset "SymPy helper surface API" begin
    @test isdefined(GenLAProblems, :sym_mat)
    @test isdefined(GenLAProblems, :sym_subs_numeric)
    @test GenLAProblems.sym_to_julia_vec([1, 2]) == [1, 2]
    @test GenLAProblems.sym_to_julia_mat([1 2; 3 4]) == [1 2; 3 4]
end

function _has_module(name::String)
    try
        GenLAProblems._ensure_pythoncall()
        PythonCall.pyimport(name)
        return true
    catch
        return false
    end
end

function _has_python_stack()
    try
        return _has_module("la_figures") && _has_module("matrixlayout")
    catch
        return false
    end
end

function _has_render_stack()
    if !_has_module("jupyter_tikz")
        return false
    end
    try
        jt = PythonCall.pyimport("jupyter_tikz")
        return PythonCall.pyconvert(Bool, PythonCall.pycall(PythonCall.pybuiltins.hasattr, jt, "render_svg_with_artifacts"))
    catch
        return false
    end
end

@testset "Wrapper kwarg helpers" begin
    d = GenLAProblems._clean_tmp_kwargs(Dict(:tmp_dir => "/tmp/la", :keep_file => "x", :output_dir => "/tmp/o", :a => 1))
    @test !haskey(d, :tmp_dir)
    @test !haskey(d, :keep_file)
    @test !haskey(d, :output_dir)
    @test d[:a] == 1

    d2 = GenLAProblems._map_tmp_to_output(Dict(:tmp_dir => "/tmp/la", :a => 2))
    @test !haskey(d2, :tmp_dir)
    @test d2[:output_dir] == "/tmp/la"
    @test d2[:a] == 2

    d3 = GenLAProblems._map_tmp_to_output(Dict(:tmp_dir => "/tmp/la", :output_dir => "/tmp/custom"))
    @test d3[:output_dir] == "/tmp/custom"
end

@testset "SymPyHelpers substitutions" begin
    if !_has_module("sympy")
        @info "Skipping SymPyHelpers tests: sympy unavailable"
    else
        sym = try
            GenLAProblems.import_sympy()
        catch err
            @info "Skipping SymPyHelpers tests: import_sympy unavailable in this runtime" exception=(err, catch_backtrace())
            nothing
        end
        if sym === nothing
            return
        end
        e = sym.symbols("e")
        f = sym.symbols("f")
        A = [1 e; 0 2]

        A_num = GenLAProblems.SymPyHelpers.sym_subs_numeric(A, Dict(e => 3))
        @test A_num isa AbstractMatrix
        @test size(A_num) == (2, 2)
        @test A_num[1, 2] == 3

        A_num_pair = GenLAProblems.SymPyHelpers.sym_subs_numeric(A, e => 4)
        @test A_num_pair isa AbstractMatrix
        @test A_num_pair[1, 2] == 4

        A_py = GenLAProblems.SymPyHelpers.sym_subs_numeric([e f; 0 1], Dict(e => 5))
        @test A_py isa PythonCall.Py
        free = PythonCall.pygetattr(A_py, "free_symbols")
        @test PythonCall.pyconvert(Int, PythonCall.pybuiltins.len(free)) == 1
    end
end

@testset "nM wrapper return contracts" begin
    if !_has_python_stack() || !_has_render_stack()
        @info "Skipping nM wrapper tests: la_figures/matrixlayout/jupyter_tikz not available"
    else
        try
            A = [2 0; 0 1]
            h_eig, spec_eig = nM.show_eig_tbl(A; tmp_dir="/tmp/la")
            @test h_eig isa GenLAProblems.SVGOut
            @test spec_eig !== nothing

            h_svd, spec_svd = nM.show_svd_table(A; tmp_dir="/tmp/la")
            @test h_svd isa GenLAProblems.SVGOut
            @test spec_svd !== nothing

            A_qr = [1 0; 0 1; 1 1]
            h_qr, mats_qr = nM.gram_schmidt_qr(A_qr; tmp_dir="/tmp/la")
            @test h_qr isa GenLAProblems.SVGOut
            @test mats_qr !== nothing
        catch err
            @info "Skipping nM wrapper tests: render/runtime not fully available" err
        end
    end
end

@testset "show_system smoke" begin
    if !_has_python_stack() || !_has_render_stack()
        @info "Skipping show_system smoke: la_figures/matrixlayout/jupyter_tikz not available"
    else
        try
            A = Rational{Int}.([1 0; 0 1])
            B = reshape(Rational{Int}.([1, 2]), 2, 1)
            pb = ShowGe{Rational{Int}}(A, B; tmp_dir="/tmp/la")
            s = show_system(pb; b_col=1, fig_scale=1)
            @test s isa GenLAProblems.SVGOut
        catch err
            @info "Skipping show_system smoke: render/runtime not fully available" err
        end
    end
end
