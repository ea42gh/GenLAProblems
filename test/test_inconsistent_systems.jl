using Test
using PythonCall
using GenLAProblems

function _py_ns_inconsistent()
    types = PythonCall.pyimport("types")
    return PythonCall.pycall(types.SimpleNamespace)
end

function _py_setattr_inconsistent(obj, name::AbstractString, value)
    PythonCall.pycall(PythonCall.pybuiltins.setattr, obj, name, value)
end

@testset "Inconsistent RHS handling" begin
    A = [1 0; 0 0]
    B = [1 0; 1 0]
    pb = ShowGE(A, B; tmp_dir="/tmp/la")
    ref!(pb; gj=true)

    @test pb.rhs_status == [:inconsistent, :consistent]
    @test pb.rhs_consistent == [false, true]
    @test pb.status == :mixed

    Xp, Xh = solutions(pb)
    @test size(Xp, 2) == 1
    @test size(Xp, 1) == size(A, 2)
    @test length(Xh) >= 0

    la = _py_ns_inconsistent()
    ml = _py_ns_inconsistent()
    seen = Ref{Dict{Symbol,Any}}(Dict{Symbol,Any}())

    _py_setattr_inconsistent(la, "backsubstitution_tex", (A, b; kwargs...) -> ["0 = 1", "\\text{No Solution}"])
    _py_setattr_inconsistent(ml, "backsubst_svg", (; kwargs...) -> (seen[] = Dict(kwargs); "<svg>cascade</svg>"))

    old_la = GenLAProblems._la_figures[]
    old_ml = GenLAProblems._matrixlayout[]
    try
        GenLAProblems._la_figures[] = la
        GenLAProblems._matrixlayout[] = ml

        h = GenLAProblems.show_backsubstitution!(pb; b_col=1)
        @test h isa GenLAProblems.SVGOut
        @test seen[][:show_cascade] == true
        @test seen[][:show_solution] == false
        @test haskey(seen[], :cascade_txt)

        sol = GenLAProblems.show_solution!(pb; b_col=1)
        @test sol == Any[]
    finally
        GenLAProblems._la_figures[] = old_la
        GenLAProblems._matrixlayout[] = old_ml
    end
end
