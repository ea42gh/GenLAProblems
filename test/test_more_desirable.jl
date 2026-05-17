using Test
using PythonCall
using GenLAProblems

function _has_module_more(name::String)
    try
        PythonCall.pyimport(name)
        return true
    catch
        return false
    end
end

function _py_setattr_more(obj, name::AbstractString, value)
    PythonCall.pycall(PythonCall.pybuiltins.setattr, obj, name, value)
end

@testset "Error-path and malformed-input contracts" begin
    if !_has_module_more("matrixlayout")
        @info "Skipping bg_for_entries error-path tests: matrixlayout unavailable"
    else
        # Malformed/unsupported bg specs should be ignored (current contract), not crash.
        @test isnothing(GenLAProblems._bg_for_entries_to_decorators("bad", [[nothing, [1 2; 3 4]]]))
        @test isnothing(GenLAProblems._bg_for_entries_to_decorators([:bad], [[nothing, [1 2; 3 4]]]))
    end

    # Non-vector specs are wrapped by normalization helper.
    one = [0, 1, [(0, 0)], "yellow!20", 1]
    @test GenLAProblems._normalize_bg_specs(one) == [one]
    @test GenLAProblems._normalize_bg_specs(nothing) == [nothing]
end

@testset "Bundle round-trip consistency under mocked backend" begin
    types = PythonCall.pyimport("types")
    la = PythonCall.pycall(types.SimpleNamespace)

    function fake_bundle(args...; kwargs...)
        spec = PythonCall.pydict(Dict("kind" => "bundle", "kws" => string(sort!(collect(keys(Dict(kwargs)))))))
        return PythonCall.pydict(Dict("spec" => spec, "svg" => "<svg>bundle</svg>"))
    end

    _py_setattr_more(la, "eig_tbl_bundle", fake_bundle)
    _py_setattr_more(la, "svd_tbl_bundle", fake_bundle)
    _py_setattr_more(la, "qr_tbl_bundle", fake_bundle)

    old_la = GenLAProblems._LAFigureSpecs[]
    try
        GenLAProblems._LAFigureSpecs[] = la

        h1, spec1 = nM.show_eig_tbl([1 0; 0 1]; tmp_dir="/tmp/la")
        raw1, spec1b = nM.eig_tbl_svg([1 0; 0 1]; tmp_dir="/tmp/la")
        raw1s = raw1 isa PythonCall.Py ? PythonCall.pyconvert(String, raw1) : String(raw1)
        @test h1 isa GenLAProblems.SVGOut
        @test raw1 isa PythonCall.Py || raw1 isa AbstractString
        @test occursin("bundle", h1.svg)
        @test occursin("bundle", raw1s)
        @test PythonCall.pyconvert(Dict, spec1) == PythonCall.pyconvert(Dict, spec1b)

        h2, spec2 = nM.show_svd_tbl([2 0; 0 1]; tmp_dir="/tmp/la")
        raw2, spec2b = nM.svd_tbl_svg([2 0; 0 1]; tmp_dir="/tmp/la")
        raw2s = raw2 isa PythonCall.Py ? PythonCall.pyconvert(String, raw2) : String(raw2)
        @test h2 isa GenLAProblems.SVGOut
        @test raw2 isa PythonCall.Py || raw2 isa AbstractString
        @test occursin("bundle", raw2s)
        @test PythonCall.pyconvert(Dict, spec2) == PythonCall.pyconvert(Dict, spec2b)
    finally
        GenLAProblems._LAFigureSpecs[] = old_la
    end
end

@testset "Type-oriented SymPy substitution behavior" begin
    if !_has_module_more("sympy")
        @info "Skipping sym_subs_numeric type tests: sympy unavailable"
    else
        sym = try
            GenLAProblems.import_sympy()
        catch err
            @info "Skipping sym_subs_numeric type tests: import_sympy unavailable in this runtime" err
            nothing
        end

        if sym !== nothing
            e = sym.symbols("e")
            A = try
                sym.Matrix(PythonCall.pylist([PythonCall.pylist([e, 1]), PythonCall.pylist([0, 2])]))
            catch err
                @info "Skipping sym_subs_numeric type tests: cannot construct SymPy matrix in this runtime" err
                nothing
            end
            if A !== nothing
                M_int = GenLAProblems.sym_subs_numeric(A, [e => 3])
                @test M_int isa AbstractMatrix
                @test M_int[1, 1] == 3

                M_float = GenLAProblems.sym_subs_numeric(A, [e => 0.5])
                @test M_float isa AbstractMatrix
                @test M_float[1, 1] isa AbstractFloat

                M_partial = GenLAProblems.sym_subs_numeric(A, Pair{Any,Any}[])
                @test M_partial isa PythonCall.Py
                free = PythonCall.pygetattr(M_partial, "free_symbols")
                @test PythonCall.pyconvert(Int, PythonCall.pybuiltins.len(free)) > 0
            end
        end
    end
end
