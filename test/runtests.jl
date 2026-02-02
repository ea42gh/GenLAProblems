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

    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    py_paths = [
        joinpath(repo, "0_ITIKZ", "la_figures"),
        joinpath(repo, "0_ITIKZ", "matrixlayout"),
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
            PythonCall.pyimport("la_figures")
            true
        catch
            @info "Skipping show_backsubstitution adjoint test: la_figures unavailable"
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
