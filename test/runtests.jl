using Test
using Random
using LinearAlgebra
using LaTeXStrings
using Aqua
using JuliaFormatter

using GenLAProblems

@testset "Package quality" begin
    Aqua.test_all(GenLAProblems; ambiguities = false, persistent_tasks = false)
end

@testset "Explicit RNG support" begin
    @test gen_qr_problem(4; rng = MersenneTwister(1234)) ==
          gen_qr_problem(4; rng = MersenneTwister(1234))
    @test gen_qr_problem(5; family = :cayley, rng = MersenneTwister(5678)) ==
          gen_qr_problem(5; family = :cayley, rng = MersenneTwister(5678))
    @test sparse_Q_matrix([2, 2]; rng = MersenneTwister(9012)) ==
          sparse_Q_matrix([2, 2]; rng = MersenneTwister(9012))
    @test gen_svd_problem(3, 2, [3, 1]; rng = MersenneTwister(3456)) ==
          gen_svd_problem(3, 2, [3, 1]; rng = MersenneTwister(3456))
end

@testset "Formatting" begin
    @test format(
        joinpath.(pkgdir(GenLAProblems), ["src", "test"]);
        overwrite = false,
        verbose = false,
    )
end

struct WrappedAlpha end
Base.show(io::IO, ::WrappedAlpha) = print(io, "\\alpha")

@testset "GenLAProblems.jl" begin
    expected_public_names = Set(
        Symbol[
            :GenLAProblems,
            :Q_2_matrix,
            :Q_3_matrix,
            :Q_4_blocks,
            :Q_4_matrix,
            :Q_matrix,
            :W_2_matrix,
            :W_3_matrix,
            :W_4_matrix,
            :W_matrix,
            :__build__,
            :__version__,
            :ca_projection_matrix,
            :e_i,
            :gen_cx_eigenproblem,
            :gen_degenerate_matrix,
            :gen_eigenproblem,
            :gen_from_jordan_form,
            :gen_full_col_rank_matrix,
            :gen_gj_matrix,
            :gen_gj_pb,
            :gen_inconsistent_gj_pb,
            :gen_inv_pb,
            :gen_ldlt_pb,
            :gen_lu_pb,
            :gen_non_diagonalizable_eigenproblem,
            :gen_particular_solution,
            :gen_permutation_matrix,
            :gen_plu_pb,
            :gen_qr_problem,
            :gen_rhs,
            :gen_svd_problem,
            :gen_symmetric_eigenproblem,
            :i_with_onecol,
            :invert_unit_lower,
            :jordan_block,
            :jordan_form,
            :lower,
            :ref_matrix,
            :rref_matrix,
            :skew_symmetric_matrix,
            :sparse_Q_matrix,
            :sparse_W_matrix,
            :symbol_vector,
            :symbols_matrix,
            :symmetric_matrix,
            :unit_lower,
        ],
    )
    @test Set(names(GenLAProblems)) == expected_public_names

    @testset "Matrix generation" begin
        @test symbol_vector(WrappedAlpha(), 1:2) ==
              [Symbol("\\alpha_1"), Symbol("\\alpha_2")]
        @test symbol_vector(L"\alpha", 1:2) == [Symbol("\\alpha_1"), Symbol("\\alpha_2")]
        @test symbols_matrix(WrappedAlpha(), 1:2, 1:3) == [
            Symbol("\\alpha_{1,1}") Symbol("\\alpha_{1,2}") Symbol("\\alpha_{1,3}")
            Symbol("\\alpha_{2,1}") Symbol("\\alpha_{2,2}") Symbol("\\alpha_{2,3}")
        ]
        @test symbols_matrix(L"\alpha", 1:2, 1:3) == [
            Symbol("\\alpha_{1,1}") Symbol("\\alpha_{1,2}") Symbol("\\alpha_{1,3}")
            Symbol("\\alpha_{2,1}") Symbol("\\alpha_{2,2}") Symbol("\\alpha_{2,3}")
        ]

        M, pivots =
            rref_matrix(4, 6, 3; maxint = 3, pivot_in_first_col = true, has_zeros = true)
        @test size(M) == (4, 6)
        @test length(pivots) == 3
        @test all(M[i, pivots[i]] == 1 for i = 1:length(pivots))

        Pvec = gen_permutation_matrix([3, 1, 2])
        @test Pvec * [10, 20, 30] == [20, 30, 10]
        @test_throws ArgumentError gen_permutation_matrix([1, 1, 2])
        @test_throws ArgumentError gen_permutation_matrix([1, 2, 3])
        Pn = gen_permutation_matrix(4)
        @test Pn * Pn' == Matrix{Int}(I, 4, 4)
        @test Pn' * Pn == Matrix{Int}(I, 4, 4)
        @test Pn != Matrix{Int}(I, 4, 4)
        @test_throws ArgumentError gen_permutation_matrix(1)

        Afull = gen_full_col_rank_matrix([2, 2], [1, 1]; maxint = 2)
        @test size(Afull) == (4, 2)
        @test rank(Afull) == 2

        pivot_cols, A =
            gen_gj_matrix(4, 6, 3; maxint = 3, pivot_in_first_col = true, has_zeros = true)
        @test size(A) == (4, 6)
        @test length(pivot_cols) == 3
        @test_throws ArgumentError rref_matrix(2, 2, 3; maxint = 2)
        @test_throws ArgumentError rref_matrix(2, 3, -1; maxint = 2)
        err_ref = try
            ref_matrix(2, 2, 3; maxint = 2)
            nothing
        catch err
            err
        end
        @test err_ref isa ArgumentError
        @test occursin(
            "ref_matrix requires r to satisfy 0 <= r <= min(m, n)",
            sprint(showerror, err_ref),
        )
        Xrhs, Brhs = gen_rhs(A, pivot_cols; maxint = 2, num_rhs = 2, has_zeros = true)
        @test size(Xrhs) == (6, 2)
        @test size(Brhs) == (4, 2)
        @test A * Xrhs == Brhs
        nonpivots = setdiff(1:size(A, 2), pivot_cols)
        @test all(Xrhs[nonpivots, :] .== 0)

        Ainc, Binc = gen_inconsistent_gj_pb(4, 6, 3; maxint = 2, num_rhs = 2)
        @test size(Ainc) == (4, 6)
        @test size(Binc) == (4, 2)
        Ar = Rational.(Ainc)
        for j in axes(Binc, 2)
            @test rank(hcat(Ar, Rational.(Binc[:, j]))) > rank(Ar)
        end
        @test_throws ArgumentError gen_inconsistent_gj_pb(3, 5, 3; maxint = 2)
        @test_throws ArgumentError gen_gj_pb(2, 2, 3; maxint = 2)
        @test_throws ArgumentError gen_lu_pb(2, 2, 3; maxint = 2)

        S, Λ, S_inv, Aeig = gen_eigenproblem([3 -1 2]; maxint = 2)
        @test size(S) == (3, 3)
        @test size(Λ) == (3, 3)
        @test size(S_inv) == (3, 3)
        @test size(Aeig) == (3, 3)

        Q, D, Asym = gen_symmetric_eigenproblem([4 1 -2]; maxint = 2)
        @test size(Q) == (3, 3)
        @test size(D) == (3, 3)
        @test size(Asym) == (3, 3)

        @test_throws ArgumentError gen_eigenproblem([1 2; 3 4]; maxint = 2)
        @test_throws ArgumentError gen_symmetric_eigenproblem([1 2; 3 4]; maxint = 2)

        @test Aeig == S * Λ * S_inv
        @test Asym == Q * D * Q'

        Sc, Λc, Sc_inv, Ac = gen_cx_eigenproblem([1 + 2im, 3]; maxint = 2)
        @test size(Λc) == (3, 3)
        @test Ac == Sc * Λc * Sc_inv
        @test Λc[1:2, 1:2] == [1 -2; 2 1]

        Smixed, Λmixed, Smixed_inv, Amixed = gen_cx_eigenproblem(
            [3, 1 // 2 + (1 // 3)im];
            maxint = 2,
            rng = MersenneTwister(7007),
        )
        @test eltype(Λmixed) == Rational{Int}
        @test Amixed == Smixed * Λmixed * Smixed_inv

        And = gen_non_diagonalizable_eigenproblem(2, 3; maxint = 2)
        @test size(And) == (3, 3)
        @test count(x -> isapprox(x, 2; atol = 1e-6, rtol = 0), eigvals(And)) >= 2

        Random.seed!(1234)
        q4_general = Q_matrix(4)
        Random.seed!(1234)
        q4_special = Q_4_matrix()
        @test q4_general == q4_special

        @test size(gen_qr_problem(4; family = :hadamard, maxint = 2)) == (4, 4)
        Random.seed!(2026)
        auto2 = gen_qr_problem(2; maxint = 2)
        Random.seed!(2026)
        pyth2 = gen_qr_problem(2; family = :pythagorean, maxint = 2)
        @test auto2 == pyth2

        @test size(gen_qr_problem(2; family = :pythagorean, maxint = 2)) == (2, 2)
        Random.seed!(2027)
        auto3 = gen_qr_problem(3; maxint = 2)
        Random.seed!(2027)
        pyth3 = gen_qr_problem(3; family = :pythagorean, maxint = 2)
        @test auto3 == pyth3
        @test size(gen_qr_problem(3; family = :pythagorean, maxint = 2)) == (3, 3)

        Random.seed!(2028)
        auto4 = gen_qr_problem(4; maxint = 2)
        Random.seed!(2028)
        pyth4 = gen_qr_problem(4; family = :pythagorean, maxint = 2)
        @test auto4 == pyth4
        @test size(gen_qr_problem(4; family = :hadamard, maxint = 2)) == (4, 4)
        @test size(gen_qr_problem(5; family = :cayley, maxint = 2)) == (5, 5)
        @test size(gen_qr_problem((2, 2); family = :sparse, maxint = 2)) == (4, 4)
        @test size(gen_qr_problem(6; maxint = 2)) == (6, 6)
        @test_throws ArgumentError gen_qr_problem(1; family = :pythagorean, maxint = 2)
        @test_throws ArgumentError gen_qr_problem((2, 2); family = :hadamard, maxint = 2)

        projA = [1 0; 0 1; 1 1]
        @test ca_projection_matrix(projA) == projA * inv(projA' * projA) * projA'
        err_proj = try
            ca_projection_matrix([1 2; 2 4])
            nothing
        catch err
            err
        end
        @test err_proj isa ArgumentError
        @test occursin(
            "ca_projection_matrix requires linearly independent columns",
            sprint(showerror, err_proj),
        )

        U1, Σ1, Vt1, Asvd1 = gen_svd_problem([2, 1], [2, 1], [3, 1, 0]; maxint = 2)
        @test U1 * Σ1 * Vt1 == Asvd1

        U2, Σ2, Vt2, Asvd2 = gen_svd_problem(
            4,
            4,
            [3, 1, 0, 0];
            left_family = :hadamard,
            right_family = :cayley,
            maxint = 2,
        )
        @test size(U2) == (4, 4)
        @test size(Vt2) == (4, 4)
        @test U2 * Σ2 * Vt2 == Asvd2

        Ainv, Ainv_inv = gen_inv_pb(4; maxint = 2)
        @test Ainv * Ainv_inv == Matrix{Int}(I, 4, 4)
        @test Ainv_inv * Ainv == Matrix{Int}(I, 4, 4)

        Ldlt, Ddlt, Adlt = gen_ldlt_pb(4; maxint = 2, rank = 3)
        @test Adlt == Ldlt * Ddlt * Ldlt'

        _, Pplu, Lplu, Uplu, Aplu = gen_plu_pb(4, 4, 3; maxint = 2)
        @test Aplu == Pplu * Lplu * Uplu
        @test Pplu != Matrix{Int}(I, 4, 4)

        Uplu0, _ = ref_matrix(4, 4, 3; maxint = 2)
        P0, L0, A0, deps0 = GenLAProblems._gen_plu_from_factors(
            Uplu0,
            3;
            maxint = 2,
            nswaps = 0,
            return_schedule = true,
        )
        @test length(deps0) == 0
        @test P0 == Matrix{Int}(I, 4, 4)
        @test A0 == P0 * L0 * Uplu0

        Uplu2, _ = ref_matrix(6, 6, 4; maxint = 2)
        P2, L2, A2, deps2 = GenLAProblems._gen_plu_from_factors(
            Uplu2,
            4;
            maxint = 2,
            nswaps = 2,
            return_schedule = true,
        )
        @test length(deps2) == 2
        @test issorted(deps2)
        @test length(unique(deps2)) == 2
        @test P2 != Matrix{Int}(I, 6, 6)
        @test A2 == P2 * L2 * Uplu2
        _, _, _, Uplu_many, Aplu_many = gen_plu_pb(6, 6, 4; maxint = 2, nswaps = 99)
        @test size(Uplu_many) == size(Aplu_many)
        @test_throws ArgumentError gen_plu_pb(2, 2, 3; maxint = 2)

        Pc, Jc, Pinvc, Ac = gen_degenerate_matrix((2, 2), (0, 1); maxint = 2)
        @test Ac == Pc * Jc * Pinvc

        Scx, Λcx, S_invcx, Acx = gen_cx_eigenproblem([-1 + 2im]; maxint = 1)
        @test Acx == Scx * Λcx * S_invcx

        @test !isdefined(GenLAProblems, :charpoly)

        Adef = Rational.(gen_non_diagonalizable_eigenproblem(2, 0; maxint = 2))
        @test rank(Adef - 2I) == 2

        J1 = jordan_block(2, 2)
        J2 = jordan_block(0, 1)
        @test_throws ArgumentError jordan_block(2, 0)
        Afrom = gen_from_jordan_form([J1, J2]; maxint = 2)
        @test_throws ArgumentError gen_degenerate_matrix((2, 0); maxint = 2)
        @test_throws ArgumentError gen_degenerate_matrix(0; maxint = 2)
    end
end
@testset "Explicit RNG direct helpers" begin
    A = [1 2; 3 4]
    @test gen_rhs(A, [1]; rng = MersenneTwister(7890)) ==
          gen_rhs(A, [1]; rng = MersenneTwister(7890))
    @test gen_particular_solution([1, 2], 3; rng = MersenneTwister(7891)) ==
          gen_particular_solution([1, 2], 3; rng = MersenneTwister(7891))
end
@testset "Explicit RNG echelon helpers" begin
    @test rref_matrix(4, 6, 3; rng = MersenneTwister(1111)) ==
          rref_matrix(4, 6, 3; rng = MersenneTwister(1111))
    @test ref_matrix(4, 6, 3; rng = MersenneTwister(2222)) ==
          ref_matrix(4, 6, 3; rng = MersenneTwister(2222))
end
@testset "Explicit RNG GE problem helpers" begin
    @test gen_gj_matrix(4, 6, 3; rng = MersenneTwister(3333)) ==
          gen_gj_matrix(4, 6, 3; rng = MersenneTwister(3333))
    @test gen_gj_pb(4, 6, 3; rng = MersenneTwister(4444)) ==
          gen_gj_pb(4, 6, 3; rng = MersenneTwister(4444))
end
@testset "Explicit RNG inconsistent GE helper" begin
    @test gen_inconsistent_gj_pb(4, 6, 3; rng = MersenneTwister(5555)) ==
          gen_inconsistent_gj_pb(4, 6, 3; rng = MersenneTwister(5555))
end

@testset "Explicit RNG factorization and basic helpers" begin
    @test gen_inv_pb(4; rng = MersenneTwister(6001)) ==
          gen_inv_pb(4; rng = MersenneTwister(6001))
    @test gen_ldlt_pb(4; rng = MersenneTwister(6002)) ==
          gen_ldlt_pb(4; rng = MersenneTwister(6002))
    @test gen_lu_pb(4, 6, 3; rng = MersenneTwister(6003)) ==
          gen_lu_pb(4, 6, 3; rng = MersenneTwister(6003))
    @test gen_plu_pb(5, 6, 3; rng = MersenneTwister(6004)) ==
          gen_plu_pb(5, 6, 3; rng = MersenneTwister(6004))
    @test gen_full_col_rank_matrix(5, 3; rng = MersenneTwister(6005)) ==
          gen_full_col_rank_matrix(5, 3; rng = MersenneTwister(6005))
    @test symmetric_matrix(4; rng = MersenneTwister(6006)) ==
          symmetric_matrix(4; rng = MersenneTwister(6006))
    @test i_with_onecol(4, 2; rng = MersenneTwister(6007)) ==
          i_with_onecol(4, 2; rng = MersenneTwister(6007))
    @test gen_permutation_matrix(4; rng = MersenneTwister(6008)) ==
          gen_permutation_matrix(4; rng = MersenneTwister(6008))
end

@testset "Explicit RNG spectral helpers" begin
    @test gen_eigenproblem([1, 2, 3]; rng = MersenneTwister(7001)) ==
          gen_eigenproblem([1, 2, 3]; rng = MersenneTwister(7001))
    @test gen_cx_eigenproblem([1 + 2im, 3]; rng = MersenneTwister(7002)) ==
          gen_cx_eigenproblem([1 + 2im, 3]; rng = MersenneTwister(7002))
    @test gen_symmetric_eigenproblem([1, 2, 3]; rng = MersenneTwister(7003)) ==
          gen_symmetric_eigenproblem([1, 2, 3]; rng = MersenneTwister(7003))
    @test gen_non_diagonalizable_eigenproblem(2, 3; rng = MersenneTwister(7004)) ==
          gen_non_diagonalizable_eigenproblem(2, 3; rng = MersenneTwister(7004))
    blocks = [jordan_block(2, 2), jordan_block(0, 1)]
    @test gen_from_jordan_form(blocks; rng = MersenneTwister(7005)) ==
          gen_from_jordan_form(blocks; rng = MersenneTwister(7005))
    @test gen_degenerate_matrix(2, (3, 1); rng = MersenneTwister(7006)) ==
          gen_degenerate_matrix(2, (3, 1); rng = MersenneTwister(7006))
end

@testset "Rank input validation" begin
    @test_throws ArgumentError rref_matrix(2.0, 3, 1)
    @test_throws ArgumentError ref_matrix(2, 3, 1.0)
    @test_throws ArgumentError gen_gj_matrix(-1, 3, 0)
    @test_throws ArgumentError gen_lu_pb(2, -3, 0)
    @test rref_matrix(3, 0, 0) == (zeros(Int64, 3, 0), Int[])
    @test rref_matrix(0, 4, 0) == (zeros(Int64, 0, 4), Int[])
    @test_throws ArgumentError rref_matrix(3, 0, 0; maxint = -1)
    @test size(ref_matrix(3, 0, 0)[1]) == (3, 0)
    @test ref_matrix(3, 4, 0; maxint = 0) == (zeros(Int64, 3, 4), Int[])
    @test size(gen_gj_matrix(3, 0, 0)[2]) == (3, 0)
    @test gen_lu_pb(3, 0, 0)[1] == Int[]
    A0, X0, B0 = gen_gj_pb(3, 0, 0)
    @test size(A0) == (3, 0)
    @test size(X0) == (0, 1)
    @test B0 == zeros(Int64, 3, 1)
end

@testset "RHS count validation" begin
    @test size(gen_rhs([1 2; 3 4], [1]; num_rhs = 0)[2]) == (2, 0)
    @test_throws ArgumentError gen_particular_solution([1], 2; num_rhs = -1)
    @test size(gen_gj_pb(2, 2, 1; num_rhs = 0)[3], 2) == 0
    @test_throws ArgumentError gen_inconsistent_gj_pb(3, 4, 1; num_rhs = 1.0)
end

@testset "maxint validation" begin
    @test_throws ArgumentError unit_lower(2; maxint = -1)
    @test_throws ArgumentError Q_matrix(2; maxint = -1)
    @test_throws ArgumentError W_3_matrix(maxint = -1)
    @test_throws ArgumentError rref_matrix(2, 2, 1; maxint = -1)
    @test_throws ArgumentError rref_matrix(2, 2, 1; maxint = 0)
    @test rref_matrix(2, 2, 1; maxint = 0, has_zeros = true)[1][1, 1] == 1
    @test_throws ArgumentError ref_matrix(2, 2, 1; maxint = 0)
    @test_throws ArgumentError gen_gj_matrix(2, 2, 1; maxint = 0)
    @test_throws ArgumentError gen_full_col_rank_matrix(3, 2; maxint = 0)
    @test_throws ArgumentError gen_inconsistent_gj_pb(3, 4, 1; maxint = 0)
    @test size(gen_gj_matrix(2, 0, 0; maxint = 0)[2]) == (2, 0)
    Aempty, Bempty = gen_inconsistent_gj_pb(2, 0, 0; maxint = 0, num_rhs = 0)
    @test size(Aempty) == (2, 0)
    @test size(Bempty) == (2, 0)
    @test_throws ArgumentError gen_particular_solution([1], 2; maxint = 1.5)
    @test unit_lower(2; maxint = 0) == Matrix{Int}(I, 2, 2)
    @test skew_symmetric_matrix(3; maxint = 0, with_zeros = true) == zeros(Int, 3, 3)
end

@testset "Dimension and index validation" begin
    @test_throws ArgumentError unit_lower(-1, 2)
    @test_throws ArgumentError e_i(0, 3)
    @test_throws ArgumentError e_i(1, 0)
    @test_throws ArgumentError i_with_onecol(3, 4)
    @test_throws ArgumentError gen_permutation_matrix(1.0)
end

@testset "Empty spectral input validation" begin
    @test_throws ArgumentError gen_eigenproblem(Int[])
    @test_throws ArgumentError gen_symmetric_eigenproblem(Int[])
    @test_throws ArgumentError gen_cx_eigenproblem(ComplexF64[])
    @test_throws ArgumentError jordan_form([])
    @test_throws ArgumentError gen_degenerate_matrix()
end

@testset "Full-column-rank input validation" begin
    @test_throws ArgumentError gen_full_col_rank_matrix(2, 3)
end

@testset "Factorization input validation" begin
    @test_throws ArgumentError gen_inv_pb(2.0)
    @test_throws ArgumentError gen_ldlt_pb(-1)
    @test_throws ArgumentError gen_ldlt_pb(3; rank = 1.5)
    @test_throws ArgumentError gen_ldlt_pb(3; maxint = 0)
    @test gen_ldlt_pb(3; maxint = 0, rank = 0)[2] == Diagonal(zeros(Int, 3))
    @test gen_inv_pb(0) == (zeros(Int, 0, 0), zeros(Int, 0, 0))
end

@testset "Block-size validation" begin
    @test size(sparse_Q_matrix(3), 1) == 3
    @test_throws ArgumentError sparse_Q_matrix((2, 0))
    @test_throws ArgumentError sparse_Q_matrix(Int[])
    @test_throws ArgumentError gen_full_col_rank_matrix((2, -1), 2)
end

@testset "QR dimension validation" begin
    @test_throws ArgumentError gen_qr_problem(0)
    @test_throws ArgumentError gen_qr_problem(-1; family = :cayley)
    @test_throws ArgumentError gen_qr_problem((2, 0); family = :sparse)
    @test_throws ArgumentError gen_qr_problem(Int[]; family = :sparse)
end

@testset "Symmetric matrix dimension validation" begin
    @test_throws ArgumentError symmetric_matrix(-1)
    @test_throws ArgumentError skew_symmetric_matrix(1.5)
    @test size(symmetric_matrix(0), 1) == 0
    @test_throws ArgumentError symmetric_matrix(2; maxint = 0)
    @test_throws ArgumentError lower(2; maxint = 0)
    @test size(lower(0; maxint = 0)) == (0, 0)
    @test size(skew_symmetric_matrix(0), 1) == 0
end

@testset "Jordan size validation" begin
    @test_throws ArgumentError jordan_block(0, 1.5)
    @test_throws ArgumentError gen_degenerate_matrix((2, 1.5))
    @test gen_degenerate_matrix(Int32(2))[2] == jordan_block(0, 2)
end

@testset "SVD input validation" begin
    @test_throws ArgumentError gen_svd_problem(0, 2, [1])
    @test_throws ArgumentError gen_svd_problem(2, (1, 0), [1])
    @test_throws ArgumentError gen_svd_problem(2, 2, Int[])
    @test_throws ArgumentError gen_svd_problem(2, 2, [1, 2, 3])
    @test_throws ArgumentError gen_svd_problem(2, 2, [-1])
    @test_throws ArgumentError gen_svd_problem(2, 2, [1 + 2im])
    @test size(gen_svd_problem(2, 3, reshape([1, 2], 1, 2))[2]) == (2, 3)
end

@testset "Complex eigenproblem input validation" begin
    @test_throws ArgumentError gen_cx_eigenproblem(1 + 2im)
    @test_throws ArgumentError gen_cx_eigenproblem(reshape([1, 2, 3, 4], 2, 2))
    @test size(gen_cx_eigenproblem(reshape([1 + 2im, 3], 1, 2))[2]) == (3, 3)
end

@testset "Jordan block collection validation" begin
    @test_throws ArgumentError jordan_form([ones(Int, 2, 3)])
    @test_throws ArgumentError jordan_form([1, 2])
    @test_throws ArgumentError jordan_form(ones(Int, 2, 2))
    @test size(jordan_form([jordan_block(0, 2)])) == (2, 2)
    mixed_jordan = jordan_form([jordan_block(1, 1), jordan_block(1 // 2, 1)])
    @test eltype(mixed_jordan) == Rational{Int}
    @test mixed_jordan == [1 0; 0 1//2]
end

@testset "RHS metadata validation" begin
    @test_throws ArgumentError gen_rhs([1 2], [0])
    @test_throws ArgumentError gen_rhs([1 2], [1, 1])
    @test_throws ArgumentError gen_rhs([1 2], [3])
    @test_throws ArgumentError gen_rhs([1, 2], [1])
    @test_throws ArgumentError gen_particular_solution([1], 0)
    @test_throws ArgumentError gen_particular_solution([1, 1], 2)
    @test gen_particular_solution(Int[], 0; num_rhs = 0) == zeros(Int, 0, 0)
    @test size(gen_rhs([1 2], [1]; num_rhs = 0, maxint = 0)[1]) == (2, 0)
    @test gen_particular_solution([1], 2; num_rhs = 0, maxint = 0) == zeros(Int, 2, 0)
    @test gen_particular_solution(Int[], 2; maxint = 0) == zeros(Int, 2, 1)
end

@testset "invert_unit_lower validation" begin
    @test invert_unit_lower([1 0; -2 1]) == [1 0; 2 1]
    @test_throws ArgumentError invert_unit_lower([1 2; 0 1])
    @test_throws ArgumentError invert_unit_lower([1 0])
    @test invert_unit_lower(zeros(Int, 0, 0)) == zeros(Int, 0, 0)
end

@testset "ca_projection_matrix validation" begin
    @test_throws ArgumentError ca_projection_matrix([1, 2])
    @test_throws ArgumentError ca_projection_matrix([1 2; 2 4])
    P = ca_projection_matrix([1 0; 0 1; 1 1])
    @test isapprox(P, P'; atol = 1e-12, rtol = 0)
    @test isapprox(P * P, P; atol = 1e-12, rtol = 0)
end
