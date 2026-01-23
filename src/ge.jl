# ------------------------------------------------------------------------------
# ------------------------------------------------------ form linear combination
# ------------------------------------------------------------------------------
raw"""
[ entries for L_show ] = form_linear_combination(s, Xh)
"""
function form_linear_combination(s, Xh)
    k    = length(s)
    expr = Vector{Any}()

    for i in 1:k
        push!(expr, s[i])
        push!(expr, Xh[:, i])
        if i < k  # Add "+" only if it's not the last term
            push!(expr, "+")
        end
    end

    return expr
end

# ==============================================================================================================

raw"""pb = ShowGe{T}(A::Matrix{T}; tmp_dir="tmp", keep_file="tmp/show\\_layout") where T <: Number
  <br>pb = ShowGe{T}(A::Matrix{T}, B::Vector{T}; tmp_dir="tmp", keep_file="tmp/show\\_layout") where T <: Number
  <br>pb = ShowGe{T}(A::Matrix{T}, B::Matrix{T}; tmp_dir="tmp", keep_file="tmp/show\\_layout") where T <: Number
  <br>pb = ShowGe{Rational{T}}(A::Matrix{T}, B::Matrix{T}; tmp_dir="tmp", keep_file="tmp/show\\_layout") where T <: Number
  <br>pb = ShowGe{Rational{T}}(A::Matrix{T}, B::Vector{T}; tmp_dir="tmp", keep_file="tmp/show\\_layout") where T <: Number
  <br>pb = ShowGe{Complex{Rational{T}}}(A::Matrix{Complex{T}}; tmp_dir="tmp", keep_file="tmp/show\\_layout") where T <: Number
  <br>pb = ShowGe{Rational{T}}(A::Matrix{T}; tmp_dir="tmp", keep_file="tmp/show\\_layout") where T <: Number
  <br>pb = ShowGe{Complex{Rational{T}}}(A::Matrix{Complex{T}}, B::Vector{Complex{T}}; tmp_dir="tmp", keep_file="tmp/show\\_layout") where T <: Number
  <br>pb = ShowGe{Complex{Rational{T}}}(A::Matrix{Complex{T}}, B::Matrix{Complex{T}}; tmp_dir="tmp", keep_file="tmp/show\\_layout") where T <: Number"""
mutable struct ShowGe{T<:Number}
    tmp_dir
    keep_file
    A
    B
    num_rhs

    matrices
    cascade
    pivot_cols
    free_cols
    desc
    pivot_list
    bg_for_entries
    ref_path_list
    basic_var
    rank
    h
    m
    xp
    xh


  function ShowGe(A::Matrix; tmp_dir="tmp", keep_file="tmp/show_layout")
	  ShowGe{eltype(A)}(A; tmp_dir=tmp_dir, keep_file=keep_file)
  end
  function ShowGe(A::Matrix, b; tmp_dir="tmp", keep_file="tmp/show_layout")
    ShowGe{eltype(A)}(A, b; tmp_dir=tmp_dir, keep_file=keep_file)
  end
  function ShowGe{T}(A::Matrix{T}; tmp_dir="tmp", keep_file="tmp/show_layout") where T <: Number
      new(tmp_dir, keep_file, A)
  end
  function ShowGe{T}(A::Matrix{T}, B::Vector{T}; tmp_dir="tmp", keep_file="tmp/show_layout") where T <: Number
      new(tmp_dir, keep_file, A,B,size(B,2))
  end
  function ShowGe{T}(A::Matrix{T}, B::Matrix{T}; tmp_dir="tmp", keep_file="tmp/show_layout") where T <: Number
      new(tmp_dir, keep_file, A,B,size(B,2))
  end

  function ShowGe{Rational{T}}(A::Matrix{T}; tmp_dir="tmp", keep_file="tmp/show_layout") where T <: Number
      new(tmp_dir, keep_file, Rational{T}.(A) )
  end
  function ShowGe{Rational{T}}(A::Matrix{T}, B::Vector{T}; tmp_dir="tmp", keep_file="tmp/show_layout") where T <: Number
      new(tmp_dir, keep_file, Rational{T}.(A),Rational{T}.(B),size(B,2))
  end
  function ShowGe{Rational{T}}(A::Matrix{T}, B::Matrix{T}; tmp_dir="tmp", keep_file="tmp/show_layout") where T <: Number
      new(tmp_dir, keep_file, Rational{T}.(A),Rational{T}.(B),size(B,2))
  end

  function ShowGe{Complex{Rational{T}}}(A::Matrix{Complex{T}}; tmp_dir="tmp", keep_file="tmp/show_layout") where T <: Number
    new(tmp_dir, keep_file, Complex{Rational{T}}.(A) )
  end
  function ShowGe{Complex{Rational{T}}}(A::Matrix{Complex{T}}, B::Vector{Complex{T}}; tmp_dir="tmp", keep_file="tmp/show_layout") where T <: Number
    new(tmp_dir, keep_file, Complex{Rational{T}}.(A),Complex{Rational{T}}.(B),size(B,2))
  end
  function ShowGe{Complex{Rational{T}}}(A::Matrix{Complex{T}}, B::Matrix{Complex{T}}; tmp_dir="tmp", keep_file="tmp/show_layout") where T <: Number
      new(tmp_dir, keep_file, Complex{Rational{T}}.(A),Complex{Rational{T}}.(B),size(B,2))
  end
end
# --------------------------------------------------------------------------------------------------------------
raw"""function ref!( pb::ShowGe{T}; N_rhs=:None, gj::Bool=false, normal\\_eq::Bool=false )  where T <: Number"""
function ref!( pb::ShowGe{T}; N_rhs=:None, gj::Bool=false, normal_eq::Bool=false )  where T <: Number
    M,N = size(pb.A)
    if isdefined( pb, :B)
       A = [pb.A pb.B]
       if N_rhs != :None
         pb.num_rhs = N_rhs
       end
    else
       A = pb.A
       pb.num_rhs = 0
    end
    if normal_eq
      pb.matrices, pb.pivot_cols, pb.desc = normal_eq_reduce_to_ref( A, n=N, gj=gj );
      sz = (N,N)
    else
      pb.matrices, pb.pivot_cols, pb.desc = reduce_to_ref( A, n=N, gj=gj );
      sz = (M,N)
    end
    pb.free_cols = filter(x -> !(x in pb.pivot_cols), 1:N)

    pb.pivot_list, pb.bg_for_entries, pb.ref_path_list, pb.basic_var = decorate_ge(pb.desc,pb.pivot_cols,sz; pivot_color="yellow!40");
    pb.rank = length( pb.pivot_cols )
    nothing
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_layout!(  pb::ShowGe{T}; array_names=nothing, show_variables=true, fig\\_scale=1 )   where T <: Number"""
function show_layout!(  pb::ShowGe{T}; array_names=nothing, show_variables=true, fig_scale=1 )   where T <: Number
    la = load_la_figures()
    rhs = isdefined(pb, :B) ? pb.B : nothing
    pb.h = la.ge_tbl_svg(pb.A, rhs;
        show_pivots=true,
        fig_scale=fig_scale,
        variable_summary=show_variables ? pb.basic_var : nothing,
        variable_colors=["red", "black"],
        array_names=array_names,
    )
    display(MIME"image/svg+xml"(), pb.h)
    return pb.h
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_system(  pb::ShowGe{T}; b_col=1, var\\_name::String="x")   where T <: Number"""
function show_system(  pb::ShowGe{T}; b_col=1, var_name::String="x")   where T <: Number
    if isdefined(pb, :B) && b_col isa Integer && 1 <= b_col <= size(pb.B, 2)
       b = pb.N[:,b_col]
    else
       b = zeros( eltype(pb.A), size(pb.A,1), 1)
    end
    tex = load_la_figures().linear_system_tex(pb.A, b, var_name=var_name)
    display(MIME"text/latex"(), tex)
    return tex
end
raw"""function show_system(  pb::ShowGe{Rational{T}}; b_col=1, var\\_name::String="x" )   where T <: Number"""
function show_system(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x" )   where T <: Number
    cnv(x) = (numerator(x),denominator(x))
    A = cnv.(pb.A)
    if isdefined(pb, :B) && b_col isa Integer && 1 <= b_col <= size(pb.B, 2)
       b = cnv.(pb.B[:,b_col])
    else
       b = cnv.(zeros( eltype(pb.A), size(A,1), 1))
    end
    tex = load_la_figures().linear_system_tex(A, b, var_name=var_name)
    display(MIME"text/latex"(), tex)
    return tex
end
raw"""function show_system(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var\\_name::String="x" )   where T <: Number"""
function show_system(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x" )   where T <: Number
    cnv(x) = (numerator(x),denominator(x))
    A = cnv.(pb.A)
    if isdefined(pb, :B) && b_col isa Integer && 1 <= b_col <= size(pb.B, 2)
       b = cnv.(pb.B[:,b_col])
    else
       b = cnv.(zeros( eltype(A), size(A,1), 1))
    end
    tex = load_la_figures().linear_system_tex(A, b, var_name=var_name)
    display(MIME"text/latex"(), tex)
    return tex
end
# --------------------------------------------------------------------------------------------------------------
raw""" cascade = create_cascade!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var\\_name::String="x" )   where T <: Number"""
function create_cascade!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x" )   where T <: Number
    pb.cascade = nothing
end
# --------------------------------------------------------------------------------------------------------------
raw""" cascade = create_cascade!(  pb::ShowGe{Rational{T}}; b_col=1, var\\_name::String="x" )   where T <: Number"""
function create_cascade!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x" )   where T <: Number
    pb.cascade = nothing
end
# --------------------------------------------------------------------------------------------------------------
raw""" cascade = create_cascade!(  pb::ShowGe{T}; b_col=nothing, var\\_name::String="x" )   where T <: Integer"""
function create_cascade!(  pb::ShowGe{T}; b_col=1, var_name::String="x" )   where T <: Integer
    pb.cascade = nothing
end
# --------------------------------------------------------------------------------------------------------------
function _encode_exact(x)
    if x isa Rational
        return (numerator(x), denominator(x))
    elseif x isa Complex{<:Rational}
        r, i = real(x), imag(x)
        return ((numerator(r), denominator(r)), (numerator(i), denominator(i)))
    end
    return x
end

function _backsub_ref(pb::ShowGe; b_col=1)
    Ab = pb.matrices[end][end]
    A = Ab[:, 1:size(pb.A, 2)]
    if isdefined(pb, :B) && b_col isa Integer && 1 <= b_col <= size(pb.B, 2)
        b = Ab[:, size(pb.A, 2) + b_col]
    else
        b = zeros(eltype(A), size(A, 1), 1)
    end
    if A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}
        A = _encode_exact.(A)
    end
    if b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}
        b = _encode_exact.(b)
    end
    return A, b
end

function _forwardsub_ref(pb::ShowGe; b_col=1)
    A = pb.A
    if isdefined(pb, :B) && b_col isa Integer && 1 <= b_col <= size(pb.B, 2)
        b = pb.B[:, b_col]
    else
        b = zeros(eltype(A), size(A, 1), 1)
    end
    if A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}
        A = _encode_exact.(A)
    end
    if b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}
        b = _encode_exact.(b)
    end
    return A, b
end

function _relabel_cascade(lines, n; var_name::String="x", param_name::String="\\alpha")
    line_list = [String(x) for x in lines]
    var_pat = Regex(string(replace(var_name, "\\" => "\\\\"), "_(\\d+)"))
    param_pat = Regex(string(replace(param_name, "\\" => "\\\\"), "_(\\d+)"))
    out = Vector{String}(undef, length(line_list))
    for (i, line) in enumerate(line_list)
        line2 = replace(line, var_pat) do m
            idx = parse(Int, m.captures[1])
            new_idx = n - idx + 1
            return string(var_name, "_", new_idx)
        end
        line2 = replace(line2, param_pat) do m
            idx = parse(Int, m.captures[1])
            new_idx = n - idx + 1
            return string(param_name, "_", new_idx)
        end
        out[i] = line2
    end
    return out
end

function _display_cascade(lines)
    tex = "\\begin{align*}\n" * join(lines, " \\\\\n") * "\n\\end{align*}"
    display(MIME"text/latex"(), tex)
    return tex
end

function _display_tex(tex)
    display(MIME"text/latex"(), tex)
    return tex
end
raw"""function show_backsubstitution!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x", fig\\_scale=1 )   where T <: Number"""
function show_backsubstitution!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Number
    A, b = _backsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A, b, var_name=var_name)
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_backsubstitution!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x", fig\\_scale=1 )   where T <: Number"""
function show_backsubstitution!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Number
    A, b = _backsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A, b, var_name=var_name)
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_backsubstitution!(  pb::ShowGe{T}; b_col=1, var_name::String="x", fig\\_scale=1 )   where T <: Integer"""
function show_backsubstitution!(  pb::ShowGe{T}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Integer
    A, b = _backsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A, b, var_name=var_name)
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_forwardsubstitution!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Number"""
function show_forwardsubstitution!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Number
    A, b = _forwardsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A[end:-1:1, end:-1:1], b[end:-1:1], var_name=var_name)
    lines = _relabel_cascade(lines, size(A, 1); var_name=var_name)
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_forwardsubstitution!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Number"""
function show_forwardsubstitution!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Number
    A, b = _forwardsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A[end:-1:1, end:-1:1], b[end:-1:1], var_name=var_name)
    lines = _relabel_cascade(lines, size(A, 1); var_name=var_name)
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_forwardsubstitution!(  pb::ShowGe{T}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Integer"""
function show_forwardsubstitution!(  pb::ShowGe{T}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Integer
    A, b = _forwardsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A[end:-1:1, end:-1:1], b[end:-1:1], var_name=var_name)
    lines = _relabel_cascade(lines, size(A, 1); var_name=var_name)
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_solution!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x", fig\\_scale=1 )   where T <: Number"""
function show_solution!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Number
    A, b = _backsub_ref(pb; b_col=b_col)
    tex = load_la_figures().standard_solution_tex(A, b, var_name=var_name)
    return _display_tex(tex)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_solution!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x", fig\\_scale=1 )   where T <: Number"""
function show_solution!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Number
    A, b = _backsub_ref(pb; b_col=b_col)
    tex = load_la_figures().standard_solution_tex(A, b, var_name=var_name)
    return _display_tex(tex)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_solution!(  pb::ShowGe{T}; b_col=1, var_name::String="x", fig\\_scale=1 )   where T <: Integer"""
function show_solution!(  pb::ShowGe{T}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Integer
    A, b = _backsub_ref(pb; b_col=b_col)
    tex = load_la_figures().standard_solution_tex(A, b, var_name=var_name)
    return _display_tex(tex)
end
# ==============================================================================================================
raw"""
    show_backsubstitution(A, b; var_name="x", fig_scale=1, tmp_dir="tmp", keep_file=nothing)

Render the back-substitution cascade for the upper-triangular system `A * x = b`
using `la_figures.backsubstitution_tex`. Works with Integer/Float as well as
exact `Rational` and `Complex{Rational}` inputs (those are converted to tuples so
SymPy reconstructs exact rationals on the Python side).
"""
function show_backsubstitution(A, b; var_name::String="x", fig_scale=1, tmp_dir="tmp", keep_file=nothing)
    A2 = (A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(A) : A
    b2 = (b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(b) : b
    lines = load_la_figures().backsubstitution_tex(A2, b2, var_name=var_name)
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
raw"""
    show_forwardsubstitution(A, b; var_name="x", fig_scale=1, tmp_dir="tmp", keep_file=nothing)

Render the forward-substitution cascade for the lower-triangular system `A * x = b`
using the la_figures backsubstitution cascade on a reversed system, then relabeling indices.
Supports Integer/Float as well as exact `Rational` and `Complex{Rational}` inputs
converted to tuples for exact SymPy reconstruction.
"""
function show_forwardsubstitution(A, b; var_name::String="x", fig_scale=1, tmp_dir="tmp", keep_file=nothing)
    A2 = (A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(A) : A
    b2 = (b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(b) : b
    lines = load_la_figures().backsubstitution_tex(A2[end:-1:1, end:-1:1], b2[end:-1:1], var_name=var_name)
    lines = _relabel_cascade(lines, size(A, 1); var_name=var_name)
    return _display_cascade(lines)
end
# ==============================================================================================================
raw"""Xp, Xh = solutions(pb::ShowGe{Complex{Rational{T}}} )   where T <: Number"""
function solutions(pb::ShowGe{Complex{Rational{T}}} )   where T <: Number
    M,N                        = size(pb.A)
    matrices, pivot_cols, desc = reduce_to_ref( pb.matrices[end][end][1:pb.rank,1:end], n = N, gj = true )

    if sum(pb.num_rhs) > 0
      Xp                         = zeros(Complex{Rational{T}}, N, sum(pb.num_rhs))
      F                          = matrices[end][end][1:pb.rank,N+1:end]
      Xp[pivot_cols,:]           = F
    else
      Xp                         = zeros(Complex{Rational{T}}, N, 1)
    end

    if length(pb.free_cols) > 0
        Xh = zeros(Complex{Rational{T}}, N, N-pb.rank)
        F  = matrices[end][end][1:pb.rank,pb.free_cols]
        for (col,row) in enumerate(pb.free_cols)  Xh[row,col] = 1  end
        Xh[pivot_cols,:] = -F
    else
        Xh = zeros(Complex{Rational{T}}, N, 1)
    end

    Xp, Xh
end
raw"""Xp, Xh = solutions(pb::ShowGe{Rational{T}} )   where T <: Number"""
function solutions(pb::ShowGe{Rational{T}} )   where T <: Number
    M,N                        = size(pb.A)
    matrices, pivot_cols, desc = reduce_to_ref( pb.matrices[end][end][1:pb.rank,1:end], n = N, gj = true )

    if sum(pb.num_rhs) > 0
      Xp                         = zeros(Rational{T}, N, sum(pb.num_rhs))
      F                          = matrices[end][end][1:pb.rank,N+1:end]
      Xp[pivot_cols,:]           = F
    else
      Xp                         = zeros(Rational{T}, N, 1)
    end

    if length(pb.free_cols) > 0
      Xh = zeros(Rational{T}, N, N-pb.rank)
      F  = matrices[end][end][1:pb.rank,pb.free_cols]
      for (col,row) in enumerate(pb.free_cols)  Xh[row,col] = 1  end
      Xh[pivot_cols,:] = -F
    else
      Xh = zeros(Rational{T}, N, 1)
    end
    Xp, Xh
end
raw"""Xp, Xh = solutions(pb::ShowGe{T} )   where T <: Number"""
function solutions(pb::ShowGe{T} )   where T <: Number
    M,N                        = size(pb.A)
    matrices, pivot_cols, desc = reduce_to_ref( pb.matrices[end][end][1:pb.rank,1:end], n = N, gj = true )

    if sum(pb.num_rhs) > 0
      Xp                         = zeros(T, N, sum(pb.num_rhs))
      F                          = matrices[end][end][1:pb.rank,N+1:end]
      Xp[pivot_cols,:]           = F
    else
      Xp                         = zeros(T, N, 1)
    end

    if length(pb.free_cols) > 0
      Xh = zeros(T, N, N-pb.rank)
      F  = matrices[end][end][1:pb.rank,pb.free_cols]
      for (col,row) in enumerate(pb.free_cols)  Xh[row,col] = 1  end
      Xh[pivot_cols,:] = -F
    else
      Xh = zeros(T, N, 1)
    end
    Xp, Xh
end
# ------------------------------------------------------------------------------------------
raw"""Xp, xH = solve!(pb::ShowGe{Complex{Rational{T}}} )   where T <: Number"""
function solve!(pb::ShowGe{Complex{Rational{T}}} )   where T <: Number
    pb.xp, pb.xh = solutions( pb )
end
raw"""solve!(pb::ShowGe{Rational{T}} )   where T <: Number"""
function solve!(pb::ShowGe{Rational{T}} )   where T <: Number
    pb.xp, pb.xh = solutions( pb )
end
raw"""Xp, Xh = solve!(pb::ShowGe{T} )   where T <: Number"""
function solve!(pb::ShowGe{T} )   where T <: Number
    pb.xp, pb.xh = solutions( pb )
end
# ==============================================================================================================
# function column_view( Xp, Xh, pivot_cols, rhs )
# end
# ==============================================================================================================
#function homogeneous_solution(pb::ShowGe{Complex{Rational{T}}}; b_col=1 )   where T <: Number)
#  N = size(pb.A,2)
#  matrices, pivot_cols, desc = reduce_to_ref( pb.matrices[end][end][:,1:N], n=N, gj=true );
#  Xh = similar(pb.A, size(pb.A,1), A - pb.rank)
#end
# ==============================================================================================================
raw"""function ge( matrices, desc, pivot_cols; Nrhs=0, formater=to_latex, pivot_list=nothing, bg_for_entries=nothing, <br>
             variable_colors=["blue","black"], pivot_colors=["blue","yellow!40"],  <br>
             ref_path_list=nothing, comment_list=[], variable_summary=nothing, array_names=nothing, <br>
             start_index=1, func=nothing, fig_scale=nothing, tmp_dir=nothing, keep_file=nothing )
"""
function julia_ge( matrices, desc, pivot_cols; Nrhs=0, formater=to_latex, pivot_list=nothing, bg_for_entries=nothing,
             variable_colors=["blue","black"], pivot_colors=["blue","yellow!40"],
             ref_path_list=nothing, comment_list=[], variable_summary=nothing, array_names=nothing,
             start_index=1, func=nothing, fig_scale=nothing, tmp_dir=nothing, keep_file=nothing )
    Ab = matrices[end][end]
    nrhs = Nrhs isa AbstractArray ? sum(Nrhs) : Nrhs
    if nrhs > 0
        A = Ab[:, 1:(end - nrhs)]
        rhs = Ab[:, (end - nrhs + 1):end]
    else
        A = Ab
        rhs = nothing
    end
    s = load_la_figures().ge_tbl_svg(A, rhs;
        fig_scale=fig_scale,
        array_names=array_names,
        variable_summary=variable_summary,
        variable_colors=variable_colors,
    )
    return s
end
"""
    SVGOut(svg::String)

Wrapper type for SVG output from `ge`, enabling rich display in IJulia.
"""
struct SVGOut
    svg::String
end

import Base: show

"""
    show(io::IO, ::MIME"image/svg+xml", x::SVGOut)

Emit the SVG payload for rich notebook display.
"""
function show(io::IO, ::MIME"image/svg+xml", x::SVGOut)
    print(io, x.svg)
end
"""
    ge(args...; kwargs...) -> SVGOut

Return an SVG wrapper for `julia_ge` suitable for notebook display.
"""
ge(args...; kwargs...) = SVGOut(julia_ge(args...; kwargs...))

# ------------------------------------------------------------------------------------------
raw"""function show_solution( matrices; var_name::String="x", tmp\\_dir=nothing )"""
function show_solution( matrices; var_name::String="x", tmp_dir=nothing )
    Ab = matrices[end][end]
    A = Ab[:, 1:(size(Ab, 2) - 1)]
    b = Ab[:, end]
    if A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}
        A = _encode_exact.(A)
    end
    if b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}
        b = _encode_exact.(b)
    end
    tex = load_la_figures().standard_solution_tex(A, b, var_name=var_name)
    return _display_tex(tex)
end
