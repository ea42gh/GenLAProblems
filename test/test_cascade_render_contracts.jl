using Test
using PythonCall
using GenLAProblems

function _py_ns_cascade()
    types = PythonCall.pyimport("types")
    return PythonCall.pycall(types.SimpleNamespace)
end

function _py_setattr_cascade(obj, name::AbstractString, value)
    PythonCall.pycall(PythonCall.pybuiltins.setattr, obj, name, value)
end

@testset "Cascade/Solution render contracts (mocked backends)" begin
    la = _py_ns_cascade()
    ml = _py_ns_cascade()

    seen_back_tex = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())
    seen_sol_tex = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())
    seen_svg = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

    function fake_backsub_tex(A, b; kwargs...)
        seen_back_tex[] = Dict(kwargs)
        return ["x_1 = 1", "x_2 = 2"]
    end

    function fake_standard_solution_tex(A, b; kwargs...)
        seen_sol_tex[] = Dict(kwargs)
        return "x = x_p + x_h"
    end

    function fake_backsubst_svg(; kwargs...)
        seen_svg[] = Dict(kwargs)
        return "<svg>cascade</svg>"
    end

    _py_setattr_cascade(la, "backsubstitution_tex", fake_backsub_tex)
    _py_setattr_cascade(la, "standard_solution_tex", fake_standard_solution_tex)
    _py_setattr_cascade(ml, "backsubst_svg", fake_backsubst_svg)

    old_la = GenLAProblems._LAFigureSpecs[]
    old_ml = GenLAProblems._matrixlayout[]

    try
        GenLAProblems._LAFigureSpecs[] = la
        GenLAProblems._matrixlayout[] = ml

        A = [1 2; 0 3]
        b = [4, 5]
        pb = ShowGE(A, b; output_dir="/tmp/la/keep")
        ref!(pb; gj=true)

        # show_backsubstitution!: should render cascade svg with expected flags.
        h1 = GenLAProblems.show_backsubstitution!(pb; b_col=1, var_name="x", fig_scale=1.4)
        @test h1 isa GenLAProblems.SVGOut
        @test occursin("cascade", h1.svg)
        @test seen_back_tex[][:var_name] == "x"
        @test seen_svg[][:show_system] == false
        @test seen_svg[][:show_cascade] == true
        @test seen_svg[][:show_solution] == false
        @test seen_svg[][:fig_scale] == 1.4
        @test seen_svg[][:output_dir] == "/tmp/la/keep"
        @test haskey(seen_svg[], :cascade_txt)

        # show_forwardsubstitution! with render_svg=false: returns plain cascade text.
        if hasmethod(show, Tuple{IO, MIME"text/latex", String})
            txt = GenLAProblems.show_forwardsubstitution!(pb; b_col=1, var_name="x", render_svg=false)
            @test txt isa AbstractString
            @test occursin("x_", txt)
        else
            @info "Skipping render_svg=false forward-substitution path: no text/latex String display method"
        end

        # show_forwardsubstitution! with render_svg=true: same svg flags as cascade path.
        h2 = GenLAProblems.show_forwardsubstitution!(pb; b_col=1, var_name="x", fig_scale=0.9, render_svg=true)
        @test h2 isa GenLAProblems.SVGOut
        @test seen_svg[][:show_system] == false
        @test seen_svg[][:show_cascade] == true
        @test seen_svg[][:show_solution] == false
        @test seen_svg[][:fig_scale] == 0.9

        # show_solution!: should route through standard_solution_tex then render solution svg.
        h3 = GenLAProblems.show_solution!(pb; b_col=1, var_name="x", fig_scale=1.1)
        @test h3 isa GenLAProblems.SVGOut
        @test seen_sol_tex[][:var_name] == "x"
        @test seen_svg[][:show_system] == false
        @test seen_svg[][:show_cascade] == false
        @test seen_svg[][:show_solution] == true
        @test seen_svg[][:fig_scale] == 1.1
        @test haskey(seen_svg[], :solution_txt)
        # Current contract: bang version does not map keep_file for solution path.
        @test !haskey(seen_svg[], :output_dir)

        h4 = GenLAProblems.show_forwardsubstitution(A, b; var_name="x", fig_scale=0.7, output_dir="/tmp/fwd")
        @test h4 isa GenLAProblems.SVGOut
        @test seen_svg[][:show_system] == false
        @test seen_svg[][:show_cascade] == true
        @test seen_svg[][:show_solution] == false
        @test seen_svg[][:fig_scale] == 0.7
        @test seen_svg[][:output_dir] == "/tmp/fwd"
    finally
        GenLAProblems._LAFigureSpecs[] = old_la
        GenLAProblems._matrixlayout[] = old_ml
    end
end

@testset "Non-bang show_solution render contract" begin
    la = _py_ns_cascade()
    ml = _py_ns_cascade()
    seen = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

    _py_setattr_cascade(la, "standard_solution_tex", (A, b; kwargs...) -> "x = x_p + x_h")
    _py_setattr_cascade(ml, "backsubst_svg", (; kwargs...) -> (seen[] = Dict(kwargs); "<svg>sol</svg>"))

    old_la = GenLAProblems._LAFigureSpecs[]
    old_ml = GenLAProblems._matrixlayout[]
    try
        GenLAProblems._LAFigureSpecs[] = la
        GenLAProblems._matrixlayout[] = ml

        mats = [[nothing, [1 0 1; 0 1 2]]]
        svg = GenLAProblems.show_solution(mats; render_svg=true, fig_scale=1.2, output_dir="/tmp/sol")
        @test svg isa GenLAProblems.SVGOut
        @test seen[][:show_solution] == true
        @test seen[][:fig_scale] == 1.2
        @test seen[][:output_dir] == "/tmp/sol"
    finally
        GenLAProblems._LAFigureSpecs[] = old_la
        GenLAProblems._matrixlayout[] = old_ml
    end
end

@testset "Non-bang show_backsubstitution render contract" begin
    la = _py_ns_cascade()
    ml = _py_ns_cascade()
    seen = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

    _py_setattr_cascade(la, "backsubstitution_tex", (A, b; kwargs...) -> ["x_2 = 1", "x_1 = 2"])
    _py_setattr_cascade(ml, "backsubst_svg", (; kwargs...) -> (seen[] = Dict(kwargs); "<svg>back</svg>"))

    old_la = GenLAProblems._LAFigureSpecs[]
    old_ml = GenLAProblems._matrixlayout[]
    try
        GenLAProblems._LAFigureSpecs[] = la
        GenLAProblems._matrixlayout[] = ml

        U = [1 2; 0 3]
        b = [4, 5]

        svg = GenLAProblems.show_backsubstitution(U, b; fig_scale=1.3, output_dir="/tmp/back")
        @test svg isa GenLAProblems.SVGOut
        @test seen[][:show_system] == false
        @test seen[][:show_cascade] == true
        @test seen[][:show_solution] == false
        @test seen[][:fig_scale] == 1.3
        @test seen[][:output_dir] == "/tmp/back"
        @test haskey(seen[], :cascade_txt)

        if hasmethod(show, Tuple{IO, MIME"text/latex", String})
            txt = GenLAProblems.show_backsubstitution(U, b; render_svg=false)
            @test txt isa AbstractString
            @test occursin("x_", txt)
        else
            @info "Skipping render_svg=false back-substitution path: no text/latex String display method"
        end
    finally
        GenLAProblems._LAFigureSpecs[] = old_la
        GenLAProblems._matrixlayout[] = old_ml
    end
end
