using Test
using PythonCall
using GenLAProblems

@testset "Additional GE helpers" begin
    @test GenLAProblems._matrices_are_strings([[nothing, ["a" "b"; "c" "d"]]])
    @test !GenLAProblems._matrices_are_strings([[nothing, [1 2; 3 4]]])

    grid = [[nothing, [1 2; 3 4]], [nothing, [5 6; 7 8]]]
    @test GenLAProblems._ge_normalize_grid(grid) == grid

    specs = [[0, 1, [[(0, 0), (1, 1)]], "yellow!20", 1]]
    @test GenLAProblems._normalize_bg_specs(specs) == specs

    one = [0, 1, [(0, 0)], "yellow!20", 1]
    @test GenLAProblems._normalize_bg_specs(one) == [one]
end

@testset "Python conversion helpers" begin
    py = GenLAProblems._ge_to_pylist(Dict("a" => Any[1, Dict("b" => 2)]))
    @test py isa PythonCall.Py
    back = PythonCall.pyconvert(Any, py)
    @test back isa AbstractDict
    @test back["a"][1] == 1
    @test back["a"][2]["b"] == 2
end

@testset "Wrapper internals" begin
    b = PythonCall.pyimport("builtins")
    d = PythonCall.pydict(Dict("spec" => "S", "svg" => "<svg/>"))
    spec, svg = GenLAProblems._bundle_result(d)
    @test PythonCall.pyconvert(String, spec) == "S"
    @test PythonCall.pyconvert(String, svg) == "<svg/>"

    types = PythonCall.pyimport("types")
    ns = PythonCall.pycall(types.SimpleNamespace)
    PythonCall.pycall(PythonCall.pybuiltins.setattr, ns, "hello", 42)
    @test PythonCall.pyconvert(Int, GenLAProblems._pygetattr_fallback(ns, :hello, "types")) == 42

    # Missing on object -> fallback module lookup.
    sn = GenLAProblems._pygetattr_fallback(ns, :SimpleNamespace, "types")
    @test PythonCall.pyconvert(Bool, PythonCall.pycall(b.callable, sn))
end
