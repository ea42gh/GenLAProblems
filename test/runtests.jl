using Test
using Random
using LinearAlgebra
using LaTeXStrings

using GenLAProblems

struct WrappedAlpha end
Base.show(io::IO, ::WrappedAlpha) = print(io, "\\alpha")

@testset "GenLAProblems.jl" begin
    @test !isdefined(GenLAProblems, :ShowGE)
    @test !isdefined(GenLAProblems, :nM)
    @test !isdefined(GenLAProblems, :PythonBridge)
    @test !isdefined(GenLAProblems, :qr_matrices_from_grid)
    @test !isdefined(GenLAProblems, :eig_matrices_from_spec)
    @test !isdefined(GenLAProblems, :svd_matrices_from_spec)
    @test !isdefined(GenLAProblems, :load_LAFigureSpecs)
    @test !isdefined(GenLAProblems, :load_matrixlayout)
    @test !isdefined(GenLAProblems, :la_version)
    @test !isdefined(GenLAProblems, :la_build)
    @test !isdefined(GenLAProblems, :ml_version)
    @test !isdefined(GenLAProblems, :ml_build)
    @test !isdefined(GenLAProblems, :materialize_python_value)
    @test !isdefined(GenLAProblems, :mm_to_px)
    @test !isdefined(GenLAProblems, :px_to_mm)
    @test !isdefined(GenLAProblems, :ensure_pythoncall!)
    @test !isdefined(GenLAProblems, :sympy)
    @test !isdefined(GenLAProblems, :_ensure_symbolics)
    @test !isdefined(GenLAProblems, :sym_to_julia_vec)
    @test !isdefined(GenLAProblems, :sym_to_julia_mat)
    @test !isdefined(GenLAProblems, :sym_subs_numeric)
    @test !isdefined(GenLAProblems, :AbstractDescription)
    @test !isdefined(GenLAProblems, :FoundPivot)
    @test !isdefined(GenLAProblems, :RequireRowExchange)
    @test !isdefined(GenLAProblems, :RequireElimination)
    @test !isdefined(GenLAProblems, :RequireScaling)
    @test !isdefined(GenLAProblems, :DoElimination)
    @test !isdefined(GenLAProblems, :DoRowExchange)
    @test !isdefined(GenLAProblems, :DoScaling)
    @test !isdefined(GenLAProblems, :Finished)
    @test !isdefined(GenLAProblems, :split_R_RHS)
    @test !isdefined(GenLAProblems, :particular_solution)
    @test !isdefined(GenLAProblems, :homogeneous_solutions)
    @test !isdefined(GenLAProblems, :find_pivot)
    @test !isdefined(GenLAProblems, :non_zero_entry)
    @test !isdefined(GenLAProblems, :interchange)
    @test !isdefined(GenLAProblems, :eliminate)
    @test !isdefined(GenLAProblems, :normal_eq_reduce_to_ref)
    @test !isdefined(GenLAProblems, :reduce_to_ref)
    @test !isdefined(GenLAProblems, :_reduce_to_ref)
    @test !isdefined(GenLAProblems, :ge_variable_type)
    @test !isdefined(GenLAProblems, :decorate_ge)
    @test !isdefined(GenLAProblems, :gram_schmidt_w)
    @test !isdefined(GenLAProblems, :normalize_columns)
    @test !isdefined(GenLAProblems, :qr_layout)
    @test !isdefined(GenLAProblems, :gram_schmidt_stable)

    @testset "Matrix generation" begin
        @test symbol_vector(WrappedAlpha(), 1:2) == [Symbol("\\alpha_1"), Symbol("\\alpha_2")]
        @test symbol_vector(L"\alpha", 1:2) == [Symbol("\\alpha_1"), Symbol("\\alpha_2")]
        @test symbols_matrix(WrappedAlpha(), 1:2, 1:3) == [
            Symbol("\\alpha_{1,1}") Symbol("\\alpha_{1,2}") Symbol("\\alpha_{1,3}");
            Symbol("\\alpha_{2,1}") Symbol("\\alpha_{2,2}") Symbol("\\alpha_{2,3}")
        ]
        @test symbols_matrix(L"\alpha", 1:2, 1:3) == [
            Symbol("\\alpha_{1,1}") Symbol("\\alpha_{1,2}") Symbol("\\alpha_{1,3}");
            Symbol("\\alpha_{2,1}") Symbol("\\alpha_{2,2}") Symbol("\\alpha_{2,3}")
        ]

        M, pivots = rref_matrix(4, 6, 3; maxint=3, pivot_in_first_col=true, has_zeros=true)
        @test size(M) == (4, 6)
        @test length(pivots) == 3
        @test all(M[i, pivots[i]] == 1 for i in 1:length(pivots))

        Pvec = gen_permutation_matrix([3, 1, 2])
        @test Pvec * [10, 20, 30] == [20, 30, 10]
        @test_throws ArgumentError gen_permutation_matrix([1, 1, 2])
        @test_throws ArgumentError gen_permutation_matrix([1, 2, 3])
        Pn = gen_permutation_matrix(4)
        @test Pn * Pn' == Matrix{Int}(I, 4, 4)
        @test Pn' * Pn == Matrix{Int}(I, 4, 4)
        @test Pn != Matrix{Int}(I, 4, 4)
        @test_throws ArgumentError gen_permutation_matrix(1)

        Afull = gen_full_col_rank_matrix([2, 2], [1, 1]; maxint=2)
        @test size(Afull) == (4, 2)
        @test rank(Afull) == 2

        pivot_cols, A = gen_gj_matrix(4, 6, 3; maxint=3, pivot_in_first_col=true, has_zeros=true)
        @test size(A) == (4, 6)
        @test length(pivot_cols) == 3
        @test_throws ArgumentError rref_matrix(2, 2, 3; maxint=2)
        @test_throws ArgumentError rref_matrix(2, 3, -1; maxint=2)
        err_ref = try
            ref_matrix(2, 2, 3; maxint=2)
            nothing
        catch err
            err
        end
        @test err_ref isa ArgumentError
        @test occursin("ref_matrix requires r to satisfy 0 <= r <= min(m, n)", sprint(showerror, err_ref))
        Xrhs, Brhs = gen_rhs(A, pivot_cols; maxint=2, num_rhs=2, has_zeros=true)
        @test size(Xrhs) == (6, 2)
        @test size(Brhs) == (4, 2)
        @test A * Xrhs == Brhs
        nonpivots = setdiff(1:size(A, 2), pivot_cols)
        @test all(Xrhs[nonpivots, :] .== 0)

        Ainc, Binc = gen_inconsistent_gj_pb(4, 6, 3; maxint=2, num_rhs=2)
        @test size(Ainc) == (4, 6)
        @test size(Binc) == (4, 2)
        Ar = Rational.(Ainc)
        for j in axes(Binc, 2)
            @test rank(hcat(Ar, Rational.(Binc[:, j]))) > rank(Ar)
        end
        @test_throws ArgumentError gen_inconsistent_gj_pb(3, 5, 3; maxint=2)
        @test_throws ArgumentError gen_gj_pb(2, 2, 3; maxint=2)
        @test_throws ArgumentError gen_lu_pb(2, 2, 3; maxint=2)

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

        projA = [1 0; 0 1; 1 1]
        @test ca_projection_matrix(projA) == projA * inv(projA' * projA) * projA'
        err_proj = try
            ca_projection_matrix([1 2; 2 4])
            nothing
        catch err
            err
        end
        @test err_proj isa ArgumentError
        @test occursin("ca_projection_matrix requires linearly independent columns", sprint(showerror, err_proj))

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
        @test Pplu != Matrix{Int}(I, 4, 4)

        Uplu0, _ = ref_matrix(4, 4, 3; maxint=2)
        P0, L0, A0, deps0 = GenLAProblems._gen_plu_from_factors(Uplu0, 3; maxint=2, nswaps=0, return_schedule=true)
        @test length(deps0) == 0
        @test P0 == Matrix{Int}(I, 4, 4)
        @test A0 == P0 * L0 * Uplu0

        Uplu2, _ = ref_matrix(6, 6, 4; maxint=2)
        P2, L2, A2, deps2 = GenLAProblems._gen_plu_from_factors(Uplu2, 4; maxint=2, nswaps=2, return_schedule=true)
        @test length(deps2) == 2
        @test issorted(deps2)
        @test length(unique(deps2)) == 2
        @test P2 != Matrix{Int}(I, 6, 6)
        @test A2 == P2 * L2 * Uplu2
        _, _, _, Uplu_many, Aplu_many = gen_plu_pb(6, 6, 4; maxint=2, nswaps=99)
        @test size(Uplu_many) == size(Aplu_many)
        @test_throws ArgumentError gen_plu_pb(2, 2, 3; maxint=2)

        Pc, Jc, Pinvc, Ac = gen_degenerate_matrix((2, 2), (0, 1); maxint=2)
        @test Ac == Pc * Jc * Pinvc

        Scx, Λcx, S_invcx, Acx = gen_cx_eigenproblem([-1 + 2im]; maxint=1)
        @test Acx == Scx * Λcx * S_invcx

        @test !isdefined(GenLAProblems, :charpoly)

        Adef = Rational.(gen_non_diagonalizable_eigenproblem(2, 0; maxint=2))
        @test rank(Adef - 2I) == 2

        J1 = jordan_block(2, 2)
        J2 = jordan_block(0, 1)
        @test_throws ArgumentError jordan_block(2, 0)
        Afrom = gen_from_jordan_form([J1, J2]; maxint=2)
        @test_throws ArgumentError gen_degenerate_matrix((2, 0); maxint=2)
        @test_throws ArgumentError gen_degenerate_matrix(0; maxint=2)
    end
end
