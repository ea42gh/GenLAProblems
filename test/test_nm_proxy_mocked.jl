using Test
using PythonCall
using GenLAProblems

function _py_ns()
    types = PythonCall.pyimport("types")
    return PythonCall.pycall(types.SimpleNamespace)
end

function _py_setattr(obj, name::AbstractString, value)
    PythonCall.pycall(PythonCall.pybuiltins.setattr, obj, name, value)
end

function _py_dict(d::Dict)
    return PythonCall.pydict(d)
end

@testset "nM wrappers with mocked backend" begin
    la = _py_ns()
    ml = _py_ns()

    seen_eig = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())
    seen_ge_tbl = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())
    seen_qr_show = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())
    seen_qr_compute = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())
    seen_render_qr = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

    function fake_bundle(args...; kwargs...)
        seen_eig[] = Dict(kwargs)
        return _py_dict(Dict(
            "spec" => _py_dict(Dict("kind" => "eig")),
            "svg" => "<svg>eig</svg>",
        ))
    end

    function fake_ge_tbl_svg(args...; kwargs...)
        seen_ge_tbl[] = Dict(kwargs)
        return "<svg>ge_tbl</svg>"
    end


    function fake_qr_svg(args...; kwargs...)
        seen_qr_show[] = Dict(kwargs)
        return "<svg>qr</svg>"
    end

    gs_obj = Any[[nothing, nothing], [nothing, nothing]]
    function fake_gs_mats(args...; kwargs...)
        seen_qr_compute[] = Dict(kwargs)
        return gs_obj
    end

    function fake_qr_spec(mats; kwargs...)
        seen_qr_compute[] = merge(seen_qr_compute[], Dict(:spec_kwargs => Dict(kwargs)))
        return _py_dict(Dict("kind" => "qr"))
    end

    function fake_render_qr_svg(; kwargs...)
        seen_render_qr[] = Dict(kwargs)
        return "<svg>qr_render</svg>"
    end

    _py_setattr(la, "eig_bundle", fake_bundle)
    _py_setattr(la, "svd_bundle", fake_bundle)
    _py_setattr(la, "qr_bundle", fake_bundle)
    _py_setattr(la, "ge_tbl_svg", fake_ge_tbl_svg)
    _py_setattr(la, "qr_svg", fake_qr_svg)
    _py_setattr(la, "gram_schmidt_qr_matrices", fake_gs_mats)
    _py_setattr(la, "qr_tbl_spec_from_matrices", fake_qr_spec)

    _py_setattr(ml, "render_qr_svg", fake_render_qr_svg)

    old_la = GenLAProblems._LAFigureSpecs[]
    old_ml = GenLAProblems._matrixlayout[]
    try
        GenLAProblems._LAFigureSpecs[] = la
        GenLAProblems._matrixlayout[] = ml

        h, spec = nM.show_eig_tbl([1 0; 0 1]; tmp_dir="/tmp/la", render_opts=Dict("crop" => "tight"))
        @test h isa GenLAProblems.SVGOut
        @test spec !== nothing
        @test haskey(seen_eig[], :output_dir)
        @test seen_eig[][:output_dir] == "/tmp/la"
        @test !haskey(seen_eig[], :tmp_dir)
        @test haskey(seen_eig[], :render_opts)

        h2 = nM.show_ge_tbl([1 0; 0 1]; tmp_dir="/tmp/la", keep_file="x", output_dir="/tmp/x", render_opts=Dict("padding" => (1, 1, 1, 1)))
        @test h2 isa GenLAProblems.SVGOut
        @test !haskey(seen_ge_tbl[], :tmp_dir)
        @test !haskey(seen_ge_tbl[], :keep_file)
        @test !haskey(seen_ge_tbl[], :output_dir)
        @test haskey(seen_ge_tbl[], :render_opts)

        h3, first = nM.show_qr([1 0; 0 1], :marker; tmp_dir="/tmp/la", render_opts=Dict("crop" => "tight"))
        @test h3 isa GenLAProblems.SVGOut
        @test first == [1 0; 0 1]
        @test haskey(seen_qr_show[], :output_dir)
        @test seen_qr_show[][:output_dir] == "/tmp/la"
        @test !haskey(seen_qr_show[], :tmp_dir)
        @test haskey(seen_qr_show[], :render_opts)

        h4, mats = nM.gram_schmidt_qr([1 0; 0 1]; tmp_dir="/tmp/la", keep_file="x", render_opts=Dict("frame" => true))
        @test h4 isa GenLAProblems.SVGOut
        @test mats isa PythonCall.Py || mats === gs_obj
        mats_j = mats isa PythonCall.Py ? PythonCall.pyconvert(Any, mats) : mats
        @test mats_j == gs_obj
        @test !haskey(seen_qr_compute[], :tmp_dir)
        @test !haskey(seen_qr_compute[], :output_dir)
        @test !haskey(seen_qr_compute[], :keep_file)
        @test haskey(seen_render_qr[], :output_dir)
        @test seen_render_qr[][:output_dir] == "/tmp/la"
        @test haskey(seen_render_qr[], :render_opts)

    finally
        GenLAProblems._LAFigureSpecs[] = old_la
        GenLAProblems._matrixlayout[] = old_ml
    end
end

@testset "nM wrappers emit deprecation guidance" begin
    la = _py_ns()
    ml = _py_ns()
    _py_setattr(la, "eig_bundle", (args...; kwargs...) -> _py_dict(Dict(
        "spec" => _py_dict(Dict("kind" => "eig")),
        "svg" => "<svg>eig</svg>",
    )))
    _py_setattr(la, "ge_tbl_svg", (args...; kwargs...) -> "<svg>ge_tbl</svg>")
    _py_setattr(la, "qr_svg", (args...; kwargs...) -> "<svg>qr</svg>")
    _py_setattr(la, "gram_schmidt_qr_matrices", (args...; kwargs...) -> Any[[nothing, nothing], [nothing, nothing]])
    _py_setattr(la, "qr_tbl_spec_from_matrices", (mats; kwargs...) -> _py_dict(Dict("kind" => "qr")))
    _py_setattr(ml, "render_qr_svg", (; kwargs...) -> "<svg>qr_render</svg>")

    old_la = GenLAProblems._LAFigureSpecs[]
    old_ml = GenLAProblems._matrixlayout[]
    try
        GenLAProblems._LAFigureSpecs[] = la
        GenLAProblems._matrixlayout[] = ml

        @test_logs (:warn, r"`nM\.show_eig_tbl` is deprecated; use `LATeachingSuite\.eig_bundle` instead\.")
            nM.show_eig_tbl([1 0; 0 1]; tmp_dir="/tmp/la")

        @test_logs (:warn, r"`nM\.show_ge_tbl` is deprecated; use `LATeachingSuite\.ge_svg \(or LATeachingSuite\.ge_bundle if you also need the spec\)` instead\.")
            nM.show_ge_tbl([1 0; 0 1]; tmp_dir="/tmp/la")

        @test_logs (:warn, r"`nM\.show_qr` is deprecated; use `LATeachingSuite\.qr_svg \(or LATeachingSuite\.qr_figure if you also need the computed matrices\)` instead\.")
            nM.show_qr([1 0; 0 1])
    finally
        GenLAProblems._LAFigureSpecs[] = old_la
        GenLAProblems._matrixlayout[] = old_ml
    end
end
