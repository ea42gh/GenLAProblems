using Test
using Random
using LinearAlgebra

using GenLAProblems

@testset "GenLAProblems.jl" begin
    @testset "Matrix generation" begin
        M, pivots = rref_matrix(4, 6, 3; maxint=3, pivot_in_first_col=true, has_zeros=true)
        @test size(M) == (4, 6)
        @test length(pivots) == 3
        @test all(M[i, pivots[i]] == 1 for i in 1:length(pivots))

        Pvec = gen_permutation_matrix([3, 1, 2])
        @test Pvec * [10, 20, 30] == [20, 30, 10]
        Pn = gen_permutation_matrix(4)
        @test Pn * Pn' == Matrix{Int}(I, 4, 4)
        @test Pn' * Pn == Matrix{Int}(I, 4, 4)

        Afull = gen_full_col_rank_matrix([2, 2], [1, 1]; maxint=2)
        @test size(Afull) == (4, 2)
        @test rank(Afull) == 2

        pivot_cols, A = gen_gj_matrix(4, 6, 3; maxint=3, pivot_in_first_col=true, has_zeros=true)
        @test size(A) == (4, 6)
        @test length(pivot_cols) == 3
        Xrhs, Brhs = gen_rhs(A, pivot_cols; maxint=2, num_rhs=2, has_zeros=true)
        @test size(Xrhs) == (6, 2)
        @test size(Brhs) == (4, 2)
        @test A * Xrhs == Brhs
        nonpivots = setdiff(1:size(A, 2), pivot_cols)
        @test all(Xrhs[nonpivots, :] .== 0)

        S, Λ, S_inv, Aeig = gen_eigenproblem([3 -1 2]; maxint=2)
        @test size(S) == (3, 3)
        @test size(Λ) == (3, 3)
        @test size(S_inv) == (3, 3)
        @test size(Aeig) == (3, 3)

        Q, D, Asym = gen_symmetric_eigenproblem([4 1 -2]; maxint=2)
        @test size(Q) == (3, 3)
        @test size(D) == (3, 3)
        @test size(Asym) == (3, 3)

        @test_throws ArgumentError gen_eigenproblem([1 2; 3 4]; maxint=2)
        @test_throws ArgumentError gen_symmetric_eigenproblem([1 2; 3 4]; maxint=2)

        Random.seed!(1234)
        q4_general = Q_matrix(4)
        Random.seed!(1234)
        q4_special = Q_4_matrix()
        @test q4_general == q4_special

        @test size(gen_qr_problem(4; family=:hadamard, maxint=2)) == (4, 4)
        Random.seed!(2026)
        auto2 = gen_qr_problem(2; maxint=2)
        Random.seed!(2026)
        pyth2 = gen_qr_problem(2; family=:pythagorean, maxint=2)
        @test auto2 == pyth2

        @test size(gen_qr_problem(2; family=:pythagorean, maxint=2)) == (2, 2)
        Random.seed!(2027)
        auto3 = gen_qr_problem(3; maxint=2)
        Random.seed!(2027)
        pyth3 = gen_qr_problem(3; family=:pythagorean, maxint=2)
        @test auto3 == pyth3
        @test size(gen_qr_problem(3; family=:pythagorean, maxint=2)) == (3, 3)

        Random.seed!(2028)
        auto4 = gen_qr_problem(4; maxint=2)
        Random.seed!(2028)
        pyth4 = gen_qr_problem(4; family=:pythagorean, maxint=2)
        @test auto4 == pyth4
        @test size(gen_qr_problem(4; family=:hadamard, maxint=2)) == (4, 4)
        @test size(gen_qr_problem(5; family=:cayley, maxint=2)) == (5, 5)
        @test size(gen_qr_problem((2, 2); family=:sparse, maxint=2)) == (4, 4)
        @test size(gen_qr_problem(6; maxint=2)) == (6, 6)
        @test_throws ArgumentError gen_qr_problem(1; family=:pythagorean, maxint=2)
        @test_throws ArgumentError gen_qr_problem((2, 2); family=:hadamard, maxint=2)

        U1, Σ1, Vt1, Asvd1 = gen_svd_problem([2, 1], [2, 1], [3, 1, 0]; maxint=2)
        @test U1 * Σ1 * Vt1 == Asvd1

        U2, Σ2, Vt2, Asvd2 = gen_svd_problem(4, 4, [3, 1, 0, 0]; left_family=:hadamard, right_family=:cayley, maxint=2)
        @test size(U2) == (4, 4)
        @test size(Vt2) == (4, 4)
        @test U2 * Σ2 * Vt2 == Asvd2

        Ainv, Ainv_inv = gen_inv_pb(4; maxint=2)
        @test Ainv * Ainv_inv == Matrix{Int}(I, 4, 4)
        @test Ainv_inv * Ainv == Matrix{Int}(I, 4, 4)

        Ldlt, Ddlt, Adlt = gen_ldlt_pb(4; maxint=2, rank=3)
        @test Adlt == Ldlt * Ddlt * Ldlt'

        _, Pplu, Lplu, Uplu, Aplu = gen_plu_pb(4, 4, 3; maxint=2)
        @test Aplu == Pplu * Lplu * Uplu

        Pc, Jc, Pinvc, Ac = gen_degenerate_matrix((2, 2), (0, 1); maxint=2)
        @test Ac == Pc * Jc * Pinvc

        Scx, Λcx, S_invcx, Acx = gen_cx_eigenproblem([-1 + 2im]; maxint=1)
        @test Acx == Scx * Λcx * S_invcx

        Adef = Rational.(gen_non_diagonalizable_eigenproblem(2, 0; maxint=2))
        @test GenLAProblems.charpoly(Adef) == GenLAProblems.charpoly([2 1 0; 0 2 0; 0 0 0])
        @test rank(Adef - 2I) == 2

        J1 = jordan_block(2, 2)
        J2 = jordan_block(0, 1)
        Afrom = gen_from_jordan_form([J1, J2]; maxint=2)
        @test GenLAProblems.charpoly(Afrom) == GenLAProblems.charpoly([2 1 0; 0 2 0; 0 0 0])
    end
end
