using Test
using PythonCall
using GenLAProblems
using LaTeXStrings

function _has_module_notebook(name::String)
    try
        PythonCall.pyimport(name)
        return true
    catch
        return false
    end
end

function _py_setattr_notebook(obj, name::AbstractString, value)
    PythonCall.pycall(PythonCall.pybuiltins.setattr, obj, name, value)
end

@testset "SymPy substitution edge-cases used in notebooks" begin
    sym = try
        GenLAProblems.import_sympy()
    catch
        nothing
    end

    if sym === nothing
        @info "Skipping SymPy substitution edge tests: sympy unavailable"
    else
        e = sym.symbols("e")
        x = sym.symbols("x")
        y = sym.symbols("y")

        # Build matrices through SymPy directly for runtime portability.
        A = sym.Matrix(PythonCall.pylist([PythonCall.pylist([1, e]), PythonCall.pylist([0, 2])]))
        B = sym.Matrix(PythonCall.pylist([PythonCall.pylist([x, 1]), PythonCall.pylist([0, y])]))

        # Common notebook pattern: Dict(symbol => value)
        A_num_float = GenLAProblems.sym_subs_numeric(A, Dict(e => 1e-3))
        @test A_num_float isa AbstractMatrix
        @test A_num_float[1, 2] isa AbstractFloat
        @test isapprox(A_num_float[1, 2], 1e-3; atol=1e-12)

        A_num_int = GenLAProblems.sym_subs_numeric(A, Dict(e => 3))
        @test A_num_int isa AbstractMatrix
        @test A_num_int[1, 2] == 3

        # Multiple substitutions in one call.
        B_num = GenLAProblems.sym_subs_numeric(B, Dict(x => 2, y => 5))
        @test B_num isa AbstractMatrix
        @test B_num[1, 1] == 2
        @test B_num[2, 2] == 5

        # Empty substitution should preserve symbolic matrix output.
        B_sym = GenLAProblems.sym_subs_numeric(B, Dict{Any,Any}())
        @test B_sym isa PythonCall.Py
    end
end

@testset "nM notebook-facing contract details" begin
    types = PythonCall.pyimport("types")
    la = PythonCall.pycall(types.SimpleNamespace)

    seen = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

    function fake_bundle(args...; kwargs...)
        seen[] = Dict(kwargs)
        return PythonCall.pydict(Dict(
            "spec" => PythonCall.pydict(Dict("kind" => "bundle-spec")),
            "svg" => "<svg>bundle</svg>",
        ))
    end

    _py_setattr_notebook(la, "svd_bundle", fake_bundle)

    old_la = GenLAProblems._LAFigureSpecs[]
    try
        GenLAProblems._LAFigureSpecs[] = la

        # Contract: show_* returns (SVGOut, spec) and tmp_dir maps to output_dir.
        h, spec = nM.show_svd_tbl([2 0; 0 1]; tmp_dir="/tmp/la", output_dir="/tmp/explicit")
        @test h isa GenLAProblems.SVGOut
        @test spec isa PythonCall.Py
        @test haskey(seen[], :output_dir)
        @test seen[][:output_dir] == "/tmp/explicit"
        @test !haskey(seen[], :tmp_dir)

        # Raw variant keeps raw svg payload in first return value.
        raw_svg, spec2 = nM.svd_tbl_svg([2 0; 0 1]; tmp_dir="/tmp/la")
        @test raw_svg isa PythonCall.Py || raw_svg isa AbstractString
        @test spec2 isa PythonCall.Py
    finally
        GenLAProblems._LAFigureSpecs[] = old_la
    end
end

@testset "GE regression guards (decorations/indexing helpers)" begin
    # Regression: Adjoint of string-like matrix cells should not attempt adjoint(::SubString).
    s = "abcdef"
    a = SubString(s, 1, 1)
    b = SubString(s, 2, 2)
    block = [a b; b a]
    adj = block'
    out = GenLAProblems._ge_block_to_list(adj)
    @test out == Any[Any[a, b], Any[b, a]]

    if try
        PythonCall.pyimport("matrixlayout")
        true
    catch
        @info "Skipping bg_for_entries regression test: matrixlayout unavailable"
        false
    end
        # Ensure range entries are preserved and passed to decorator conversion.
        bg_specs = [[0, 1, [[(0, 0), (2, 3)], [(1, 1), (2, 2)]], "yellow!20", 1]]
        mats_raw = [[nothing, [1 2 3 4; 5 6 7 8; 9 10 11 12]]]
        decs = GenLAProblems._bg_for_entries_to_decorators(bg_specs, mats_raw)
        @test decs !== nothing
        @test length(decs) == 1
        @test decs[1]["grid"] == (0, 1)
    end

@testset "QR matrices from grid are Julia matrices for l_show" begin
    if !_has_module_notebook("jupyter_tikz") || Sys.which("latexmk") === nothing
        @info "Skipping QR matrices from grid regression test: jupyter_tikz/latexmk unavailable"
    else
        GenLAProblems.ensure_pythoncall!()
        A = [1 -2 -1 1; -1 2 3 -1; 1 0 1 -3; -1 0 1 -1]
        h, mats = nM.gram_schmidt_qr(A; fig_scale=1.2)
        qr = qr_matrices_from_grid(mats)

        @test qr.Q isa AbstractMatrix
        @test qr.R isa AbstractMatrix
        @test !(qr.Q isa PythonCall.Py)
        @test !(qr.R isa PythonCall.Py)

        latex = l_show(L"A = Q R = ", qr.Q, qr.R)
        @test occursin("Q", String(latex))
        @test occursin("R", String(latex))
    end
end
end
