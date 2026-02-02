using Test
using PythonCall
using GenLAProblems

function _has_module_rt(name::String)
    try
        PythonCall.pyimport(name)
        return true
    catch
        return false
    end
end

function _py_setattr_rt(obj, name::AbstractString, value)
    PythonCall.pycall(PythonCall.pybuiltins.setattr, obj, name, value)
end

@testset "Module cache helpers" begin
    types = PythonCall.pyimport("types")
    la = PythonCall.pycall(types.SimpleNamespace)
    ml = PythonCall.pycall(types.SimpleNamespace)
    old_la = GenLAProblems._la_figures[]
    old_ml = GenLAProblems._matrixlayout[]
    try
        GenLAProblems._la_figures[] = la
        GenLAProblems._matrixlayout[] = ml
        @test GenLAProblems.load_la_figures() === la
        @test GenLAProblems.load_matrixlayout() === ml
        @test GenLAProblems.load_la_figures() === la
        @test GenLAProblems.load_matrixlayout() === ml
    finally
        GenLAProblems._la_figures[] = old_la
        GenLAProblems._matrixlayout[] = old_ml
    end
end

@testset "SVG helper negative contract" begin
    if !_has_module_rt("IPython.display")
        @info "Skipping SVG helper negative tests: IPython.display unavailable"
    else
        err = try
            GenLAProblems.py_show_svg(123)
            nothing
        catch e
            e
        end
        @test err !== nothing
        @test occursin("py_show_svg expects", sprint(showerror, err))
    end
end

@testset "show_solution render_svg=false path" begin
    if !hasmethod(show, Tuple{IO, MIME"text/latex", String})
        @info "Skipping show_solution render_svg=false test: no text/latex String display method in this runtime"
    else
    types = PythonCall.pyimport("types")
    la = PythonCall.pycall(types.SimpleNamespace)
    old_la = GenLAProblems._la_figures[]
    try
        _py_setattr_rt(la, "standard_solution_tex", (A, b; kwargs...) -> "x_1 = 1")
        GenLAProblems._la_figures[] = la
        mats = [[nothing, [1 0 1; 0 1 2]]]
        tex = GenLAProblems.show_solution(mats; render_svg=false)
        @test tex isa AbstractString
        @test occursin("x_1 = 1", tex)
    finally
        GenLAProblems._la_figures[] = old_la
    end
    end
end
