using Test
using PythonCall

# Ensure Python modules in local workspace are importable during tests.
let
    # Pin PythonCall to the same interpreter used for local notebook/runtime checks.
    if !haskey(ENV, "JULIA_PYTHONCALL_EXE")
        ENV["JULIA_PYTHONCALL_EXE"] = get(ENV, "PYTHON", Sys.which("python3"))
    end
    # Avoid CondaPkg creating an isolated Python during tests.
    ENV["JULIA_CONDAPKG_BACKEND"] = get(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
    ENV["CONDAPKG_BACKEND"] = get(ENV, "CONDAPKG_BACKEND", "Null")

    repo = normpath(joinpath(@__DIR__, "..", ".."))
    py_paths = [
        joinpath(repo, ".pydeps"),
        joinpath(repo, "LAFigureSpecs"),
        joinpath(repo, "matrixlayout"),
    ]
    existing = get(ENV, "PYTHONPATH", "")
    parts = isempty(existing) ? String[] : split(existing, ':')
    for p in reverse(py_paths)
        if !(p in parts)
            pushfirst!(parts, p)
        end
    end
    ENV["PYTHONPATH"] = join(parts, ':')
end

using GenLAProblems

@testset "GenLAProblems.jl" begin
    @testset "Matrix generation" begin
        M, pivots = rref_matrix(4, 6, 3; maxint=3, pivot_in_first_col=true, has_zeros=true)
        @test size(M) == (4, 6)
        @test length(pivots) == 3
        @test all(M[i, pivots[i]] == 1 for i in 1:length(pivots))

        pivot_cols, A = gen_gj_matrix(4, 6, 3; maxint=3, pivot_in_first_col=true, has_zeros=true)
        @test size(A) == (4, 6)
        @test length(pivot_cols) == 3

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

        @test size(gen_qr_problem_hadamard(4; maxint=2)) == (4, 4)
        @test size(gen_qr_problem(2; family=:pythagorean, maxint=2)) == (2, 2)
        @test size(gen_qr_problem(3; family=:pythagorean, maxint=2)) == (3, 3)
        @test size(gen_qr_problem(4; family=:hadamard, maxint=2)) == (4, 4)
        @test size(gen_qr_problem(5; family=:cayley, maxint=2)) == (5, 5)
        @test size(gen_qr_problem((2, 2); family=:sparse, maxint=2)) == (4, 4)
        @test size(gen_qr_problem(3; maxint=2)) == (3, 3)
        @test size(gen_qr_problem(6; maxint=2)) == (6, 6)
        @test_throws ArgumentError gen_qr_problem(1; family=:pythagorean, maxint=2)
        @test_throws ArgumentError gen_qr_problem((2, 2); family=:hadamard, maxint=2)

        U1, Σ1, Vt1, Asvd1 = gen_svd_problem([2, 1], [2, 1], [3, 1, 0]; maxint=2)
        @test U1 * Σ1 * Vt1 == Asvd1

        U2, Σ2, Vt2, Asvd2 = gen_svd_problem(4, 4, [3, 1, 0, 0]; left_family=:hadamard, right_family=:cayley, maxint=2)
        @test size(U2) == (4, 4)
        @test size(Vt2) == (4, 4)
        @test U2 * Σ2 * Vt2 == Asvd2
    end

    @testset "Reduction helpers" begin
        A = [1 2 3; 2 4 6; 1 1 1]
        matrices, pivots, _ = reduce_to_ref(A; gj=true)
        @test size(matrices[end][end], 2) == size(A, 2)
        @test length(pivots) == 2
    end

    @testset "Adjoint inputs" begin
        A = Rational.( [1 2; 3 4] )
        Aadj = A'
        matrices, _, _ = reduce_to_ref(Aadj; gj=true)
        @test size(matrices[end][end]) == (2, 2)

        Q, R = gram_schmidt_stable(Aadj)
        @test size(Q) == (2, 2)
        @test size(R) == (2, 2)

        @test charpoly(Aadj) == charpoly(Matrix(Aadj))

        matrices2, _, _ = normal_eq_reduce_to_ref(Aadj; gj=false)
        @test size(matrices2[end][end]) == (2, 2)

        W = gram_schmidt_w(Aadj)
        @test size(W) == (2, 2)
        @test eltype(W) <: Rational

        layout = qr_layout(Aadj)
        @test layout !== nothing

        P = ca_projection_matrix(Aadj)
        @test size(P) == (2, 2)
        @test P == P'

        U = [Rational(1) Rational(2); 0 Rational(3)]
        b = Rational.( [5, 6] )
        if try
            GenLAProblems._ensure_pythoncall()
            PythonCall.pyimport("LAFigureSpecs")
            PythonCall.pyimport("jupyter_tikz")
            Sys.which("latexmk") !== nothing
        catch
            @info "Skipping show_backsubstitution adjoint test: LAFigureSpecs/jupyter_tikz/latexmk unavailable"
            false
        end
            @test show_backsubstitution(U', b) !== nothing
        end

        L = [Rational(1) 0; Rational(2) Rational(3)]
        if hasmethod(show, Tuple{IO, MIME"text/latex", String})
            @test show_forwardsubstitution(L', b; render_svg=false) !== nothing
        else
            @info "Skipping show_forwardsubstitution render_svg=false adjoint test: no text/latex String display method"
        end
    end

    @testset "Solve helpers" begin
        R_RHS = [1 0 2 5; 0 1 -1 3]
        R, RHS = split_R_RHS(R_RHS, 1)
        @test R == [1 0 2; 0 1 -1]
        @test RHS == [5; 3;;]

        R2, groups = split_R_RHS([1 2 10 11 20; 3 4 12 13 21], [2, 1])
        @test R2 == [1 2; 3 4]
        @test length(groups) == 2
        @test groups[1] == [10 11; 12 13]
        @test groups[2] == [20; 21;;]
        @test_throws ArgumentError split_R_RHS(R_RHS, -1)
        @test_throws ArgumentError split_R_RHS(R_RHS, [1, -1])
        @test_throws ArgumentError split_R_RHS(R_RHS, [10, 10])

        pivots = [1, 2]
        Xp = particular_solution(R, RHS, pivots)
        @test Xp == [5; 3; 0;;]

        H = homogeneous_solutions(R, pivots)
        @test H == [-2; 1; 1;;]

        A = [0 2 0; 0 0 3; 4 5 6]
        @test GenLAProblems.find_pivot(A, 1, 1) == 3
        @test GenLAProblems.find_pivot(A, 2, 2) == 3
        @test GenLAProblems.find_pivot(A, 2, 1) == 3
        @test GenLAProblems.find_pivot(A, 3, 2) == 3

        B = [1 0 0; 0 0 2; 0 3 0]
        @test GenLAProblems.non_zero_entry(B, 1, 2, false)
        @test GenLAProblems.non_zero_entry(B, 2, 2, false)
        @test GenLAProblems.non_zero_entry(B, 2, 2, true)

        C = [1 2; 3 4]
        GenLAProblems.interchange(C, 1, 2)
        @test C == [3 4; 1 2]
        GenLAProblems.eliminate(C, 1, 2, -2)
        @test C == [3 4; -5 -6]
    end

    @testset "Normal equation layout uses computed matrices" begin
        if try
            GenLAProblems._ensure_pythoncall()
            PythonCall.pyimport("LAFigureSpecs")
            PythonCall.pyimport("matrixlayout")
            PythonCall.pyimport("jupyter_tikz")
            Sys.which("latexmk") !== nothing
        catch
            @info "Skipping normal_eq layout test: Python deps/render stack/latexmk unavailable"
            false
        end
            A = [1 2 3; 4 5 6]
            b = [1, 1]
            pb = ShowGE{Rational{Int}}(A, b)
            ref!(pb; normal_eq=true)
            # Rendering should use the normal-equation matrix stack.
            h = show_layout!(pb; array_names=["A", ["Aᵀ", "AᵀA"]], fig_scale=1, render_opts=Dict("crop"=>"tight"))
            @test h !== nothing
        end
    end
end

include("test_python_wrappers.jl")
include("test_ge_helpers.jl")
include("test_bg_decorators.jl")
include("test_nm_proxy_mocked.jl")
include("test_wrapper_and_ge_more.jl")
include("test_additional_coverage.jl")
include("test_more_desirable.jl")
include("test_runtime_surface.jl")
include("test_notebook_contracts_and_regressions.jl")
include("test_cascade_render_contracts.jl")
include("test_show_system_contracts.jl")
include("test_sym_ops_and_extractors.jl")
include("test_inconsistent_systems.jl")
