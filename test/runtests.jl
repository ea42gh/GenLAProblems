using GenLAProblems
using Test

@testset "GenLAProblems.jl" begin
    @testset "Matrix generation" begin
        M, pivots = rref_matrix(4, 6, 3; maxint=3, pivot_in_first_col=true, has_zeros=true)
        @test size(M) == (4, 6)
        @test length(pivots) == 3
        @test all(M[i, pivots[i]] == 1 for i in 1:length(pivots))

        pivot_cols, A = gen_gj_matrix(4, 6, 3; maxint=3, pivot_in_first_col=true, has_zeros=true)
        @test size(A) == (4, 6)
        @test length(pivot_cols) == 3
    end

    @testset "Reduction helpers" begin
        A = [1 2 3; 2 4 6; 1 1 1]
        matrices, pivots, _ = reduce_to_ref(A; gj=true)
        @test size(matrices[end][end], 2) == size(A, 2)
        @test length(pivots) == 2
    end

    @testset "Solve helpers" begin
        R_RHS = [1 0 2 5; 0 1 -1 3]
        R, RHS = split_R_RHS(R_RHS, 1)
        @test R == [1 0 2; 0 1 -1]
        @test RHS == [5; 3]

        pivots = [1, 2]
        Xp = particular_solution(R, RHS, pivots)
        @test Xp == [5; 3; 0]

        H = homogeneous_solutions(R, pivots)
        @test H == [-2; 1; 1]

        A = [0 2 0; 0 0 3; 4 5 6]
        @test find_pivot(A, 1, 1) == 3
        @test find_pivot(A, 2, 2) == 3
        @test find_pivot(A, 2, 1) == 3
        @test find_pivot(A, 3, 2) == -1
        @test find_diag_pivot(A, 1, 1) == 3
        @test find_diag_pivot(A, 2, 2) == 3
        @test find_diag_pivot(A, 3, 3) == 3
        @test find_diag_pivot(zeros(2, 2), 1, 1) == -1

        B = [1 0 0; 0 0 2; 0 3 0]
        @test non_zero_entry(B, 1, 2, false)
        @test !non_zero_entry(B, 2, 2, false)
        @test non_zero_entry(B, 2, 2, true)

        C = [1 2; 3 4]
        interchange(C, 1, 2)
        @test C == [3 4; 1 2]
        eliminate(C, 1, 2, -2)
        @test C == [3 4; -5 -6]
    end
end
