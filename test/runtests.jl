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
end
