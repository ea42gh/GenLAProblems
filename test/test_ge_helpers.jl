using Test
using PythonCall
using GenLAProblems

@testset "GE helper normalization" begin
    mats = [[1 2; 3 4], [5 6; 7 8]]
    normed = GenLAProblems._ge_normalize_grid(mats)
    @test normed == [[[1 2; 3 4]], [[5 6; 7 8]]]

    grid = [[nothing, [1 2; 3 4]], [:none, [5 6; 7 8]]]
    lists = GenLAProblems._ge_grid_to_lists(grid)
    @test lists[1][1] === nothing
    @test lists[2][1] === nothing
    @test lists[1][2] == Any[Any[1, 2], Any[3, 4]]
    @test lists[2][2] == Any[Any[5, 6], Any[7, 8]]
end

@testset "GE helper spec normalization and payload conversion" begin
    GenLAProblems.ensure_pythoncall!()

    one_spec = [0, 1, [[(0, 0), (1, 1)]], "yellow!20", 1]
    @test GenLAProblems._normalize_bg_specs(one_spec) == [one_spec]
    @test GenLAProblems._normalize_bg_specs([one_spec]) == [one_spec]

    mats_raw = [[nothing, [1 2; 3 4]]]
    @test GenLAProblems._matrix_shape_from_grid(mats_raw, 0, 1) == (2, 2)
    @test GenLAProblems._matrix_shape_from_grid(mats_raw, 0, 0) == (0, 0)

    payload = GenLAProblems._ge_pyify_payload(
        [[nothing, Any[Any[1, 2], Any[3, 4]]]],
        nothing, nothing, nothing, [], nothing, nothing, nothing, nothing, nothing
    )
    @test payload.mats isa PythonCall.Py
    @test payload.comment_list isa PythonCall.Py
end

@testset "GE helper block conversion" begin
    @test GenLAProblems._ge_block_to_list(nothing) === nothing
    @test GenLAProblems._ge_block_to_list(:none) === nothing
    @test GenLAProblems._ge_block_to_list("x") == "x"
    @test GenLAProblems._ge_block_to_list([1 2; 3 4]) == Any[Any[1, 2], Any[3, 4]]
end

@testset "n_rhs divider vector compatibility" begin
    A = Rational.( [1 0; 0 1] )
    B = Rational.( [2 3; 4 5] )
    pb = ShowGE(A, B)
    ref!(pb; gj=false, n_rhs=[1, 1])

    @test pb.n_rhs == [1, 1]
    pb2 = ShowGE(A, B; output_dir="/tmp/alt")
    @test pb2.tmp_dir == "/tmp/alt"
    @test length(pb.rhs_status) == size(B, 2)
    @test length(pb.rhs_consistent) == size(B, 2)
    @test all(s -> s == :consistent, pb.rhs_status)
end
