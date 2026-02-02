using Test
using PythonCall
using GenLAProblems

function _has_module_bg(name::String)
    try
        PythonCall.pyimport(name)
        return true
    catch
        return false
    end
end

@testset "_bg_for_entries_to_decorators semantics" begin
    if !_has_module_bg("matrixlayout")
        @info "Skipping bg_for_entries decorator tests: matrixlayout unavailable"
    else
        mats = [
            [nothing, [1 2; 3 4]],
            [[1 0; 0 1], [5 6; 7 8]],
        ]

        # Mixed selectors: box range + single entry.
        specs = [
            [0, 1, [[(0, 0), (1, 1)]], "yellow!35", 1],
            [1, 1, [(0, 1)], "blue!20", 1],
        ]
        decs = GenLAProblems._bg_for_entries_to_decorators(specs, mats)
        @test decs isa Vector
        @test length(decs) == 2
        @test decs[1]["grid"] == (0, 1)
        @test decs[2]["grid"] == (1, 1)
        @test length(decs[1]["entries"]) == 1
        @test length(decs[2]["entries"]) == 1

        callable = PythonCall.pybuiltins.callable
        @test PythonCall.pyconvert(Bool, PythonCall.pycall(callable, decs[1]["decorator"]))
        @test PythonCall.pyconvert(Bool, PythonCall.pycall(callable, decs[2]["decorator"]))

        # Single spec should normalize to a one-item decorators list.
        one = [0, 1, [[(0, 0), (0, 0)]], "red!15", 1]
        decs_one = GenLAProblems._bg_for_entries_to_decorators(one, mats)
        @test decs_one isa Vector
        @test length(decs_one) == 1
        @test decs_one[1]["grid"] == (0, 1)

        # Nothing stays nothing.
        @test isnothing(GenLAProblems._bg_for_entries_to_decorators(nothing, mats))
    end
end
