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
    ge_tbl_svg = _pygetattr(la, :ge_tbl_svg)
    pb.h = _pycall(ge_tbl_svg, pb.A, rhs;
        show_pivots=true,
        fig_scale=fig_scale,
        variable_summary=show_variables ? pb.basic_var : nothing,
        variable_colors=["red", "black"],
        array_names=array_names,
    )
    _ensure_pythoncall()
    svg_str = Base.invokelatest(PythonCall.pyconvert, String, pb.h)
    svg = SVGOut(svg_str)
    pb.h = svg
    return svg
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_system(  pb::ShowGe{T}; b_col=1, var\\_name::String="x")   where T <: Number"""
function show_system(  pb::ShowGe{T}; b_col=1, var_name::String="x")   where T <: Number
    if isdefined(pb, :B) && b_col isa Integer && 1 <= b_col <= size(pb.B, 2)
       b = pb.B[:,b_col]
    else
       b = zeros( eltype(pb.A), size(pb.A,1), 1)
    end
    la = load_la_figures()
    linear_system_tex = _pygetattr(la, :linear_system_tex)
    tex = _pycall(linear_system_tex, pb.A, b; var_name=var_name)
    _ensure_pythoncall()
    tex = Base.invokelatest(PythonCall.pyconvert, String, tex)
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
    la = load_la_figures()
    linear_system_tex = _pygetattr(la, :linear_system_tex)
    tex = _pycall(linear_system_tex, A, b; var_name=var_name)
    _ensure_pythoncall()
    tex = Base.invokelatest(PythonCall.pyconvert, String, tex)
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
    la = load_la_figures()
    linear_system_tex = _pygetattr(la, :linear_system_tex)
    tex = _pycall(linear_system_tex, A, b; var_name=var_name)
    _ensure_pythoncall()
    tex = Base.invokelatest(PythonCall.pyconvert, String, tex)
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

function _rational_str(r::Rational)
    return string(numerator(r), "/", denominator(r))
end

function _complex_rational_str(z::Complex{<:Rational})
    re = real(z)
    im = imag(z)
    if im == 0
        return _rational_str(re)
    end
    im_str = _rational_str(abs(im))
    if re == 0
        return string(im < 0 ? "-" : "", im_str, "*I")
    end
    sign = im < 0 ? "-" : "+"
    return string(_rational_str(re), " ", sign, " ", im_str, "*I")
end

function _encode_exact_vector(b::AbstractVector)
    out = Vector{Any}(undef, length(b))
    for i in eachindex(b)
        val = b[i]
        if val isa Rational
            out[i] = _rational_str(val)
        elseif val isa Complex{<:Rational}
            out[i] = _complex_rational_str(val)
        else
            out[i] = val
        end
    end
    return out
end

function _rhs_vector(b, b_col)
    if b isa AbstractMatrix
        if size(b, 2) == 1
            return vec(b)
        end
        if b_col isa Integer && 1 <= b_col <= size(b, 2)
            return b[:, b_col]
        end
        return b[:, 1]
    end
    return b
end

function _backsub_ref(pb::ShowGe; b_col=1)
    Ab = pb.matrices[end][end]
    if Ab isa AbstractArray{<:AbstractString} || any(x -> x isa AbstractString, Ab)
        gj = false
        if isdefined(pb, :desc)
            for d in pb.desc
                if hasproperty(d, :gj) && getproperty(d, :gj) === true
                    gj = true
                    break
                end
            end
        end
        if isdefined(pb, :B)
            Ab_full = [pb.A pb.B]
        else
            Ab_full = pb.A
        end
        mats, _, _ = reduce_to_ref(Ab_full, n=size(pb.A, 2), gj=gj)
        Ab = mats[end][end]
    end
    A = Ab[:, 1:size(pb.A, 2)]
    if isdefined(pb, :B) && b_col isa Integer && 1 <= b_col <= size(pb.B, 2)
        b = Ab[:, size(pb.A, 2) + b_col]
    else
        b = zeros(eltype(A), size(A, 1), 1)
    end
    b = _rhs_vector(b, b_col)
    if A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}
        A = _encode_exact.(A)
    end
    if b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}
        if b isa AbstractVector
            b = _encode_exact_vector(b)
        else
            b = _encode_exact.(b)
        end
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
    b = _rhs_vector(b, b_col)
    if A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}
        A = _encode_exact.(A)
    end
    if b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}
        if b isa AbstractVector
            b = _encode_exact_vector(b)
        else
            b = _encode_exact.(b)
        end
    end
    return A, b
end

function _relabel_cascade(lines, n; var_name::String="x", param_name::String="\\alpha")
    if !(lines isa AbstractVector{<:AbstractString})
        _ensure_pythoncall()
        lines = Base.invokelatest(PythonCall.pyconvert, Vector{String}, lines)
    end
    line_list = [String(x) for x in lines]
    var_pat = Regex(string(replace(var_name, "\\" => "\\\\"), "_(\\d+)"))
    param_pat = Regex(string(replace(param_name, "\\" => "\\\\"), "_(\\d+)"))
    out = Vector{String}(undef, length(line_list))
    for (i, line) in enumerate(line_list)
        line2 = replace(line, var_pat => (s -> begin
            m = match(var_pat, s)
            idx = parse(Int, m.captures[1])
            new_idx = n - idx + 1
            string(var_name, "_", new_idx)
        end))
        line2 = replace(line2, param_pat => (s -> begin
            m = match(param_pat, s)
            idx = parse(Int, m.captures[1])
            new_idx = n - idx + 1
            string(param_name, "_", new_idx)
        end))
        out[i] = line2
    end
    return out
end

function _display_cascade(lines)
    if !(lines isa AbstractVector{<:AbstractString})
        _ensure_pythoncall()
        lines = PythonCall.pyconvert(Vector{String}, lines)
    end
    tex = join(lines, "\n")
    display(MIME"text/latex"(), tex)
    return tex
end

function _render_backsubst_svg(lines; fig_scale=nothing, tmp_dir=nothing, keep_file=nothing)
    ml = load_matrixlayout()
    backsubst_svg = _pygetattr(ml, :backsubst_svg)
    kwargs = Dict{Symbol, Any}()
    kwargs[:cascade_txt] = _ge_to_pylist(lines)
    kwargs[:show_system] = false
    kwargs[:show_cascade] = true
    kwargs[:show_solution] = false
    if fig_scale !== nothing
        kwargs[:fig_scale] = fig_scale
    end
    if tmp_dir !== nothing
        kwargs[:tmp_dir] = tmp_dir
    end
    svg = _pycall(backsubst_svg; kwargs...)
    return _show_svg(svg)
end

function _render_solution_svg(solution_tex; fig_scale=nothing, tmp_dir=nothing)
    ml = load_matrixlayout()
    backsubst_svg = _pygetattr(ml, :backsubst_svg)
    kwargs = Dict{Symbol, Any}()
    kwargs[:solution_txt] = solution_tex
    kwargs[:show_system] = false
    kwargs[:show_cascade] = false
    kwargs[:show_solution] = true
    if fig_scale !== nothing
        kwargs[:fig_scale] = fig_scale
    end
    if tmp_dir !== nothing
        kwargs[:tmp_dir] = tmp_dir
    end
    svg = _pycall(backsubst_svg; kwargs...)
    return _show_svg(svg)
end

function _display_tex(tex)
    if !(tex isa AbstractString)
        _ensure_pythoncall()
        tex = Base.invokelatest(PythonCall.pyconvert, String, tex)
    end
    display(MIME"text/latex"(), tex)
    return tex
end
raw"""function show_backsubstitution!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x", fig\\_scale=1 )   where T <: Number"""
function show_backsubstitution!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Number
    A, b = _backsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A, b, var_name=var_name)
    return _render_backsubst_svg(lines; fig_scale=fig_scale, tmp_dir=pb.tmp_dir, keep_file=pb.keep_file)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_backsubstitution!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x", fig\\_scale=1 )   where T <: Number"""
function show_backsubstitution!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Number
    A, b = _backsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A, b, var_name=var_name)
    return _render_backsubst_svg(lines; fig_scale=fig_scale, tmp_dir=pb.tmp_dir, keep_file=pb.keep_file)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_backsubstitution!(  pb::ShowGe{T}; b_col=1, var_name::String="x", fig\\_scale=1 )   where T <: Integer"""
function show_backsubstitution!(  pb::ShowGe{T}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Integer
    A, b = _backsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A, b, var_name=var_name)
    return _render_backsubst_svg(lines; fig_scale=fig_scale, tmp_dir=pb.tmp_dir, keep_file=pb.keep_file)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_forwardsubstitution!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x", fig_scale=1, render_svg=true )   where T <: Number"""
function show_forwardsubstitution!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x", fig_scale=1, render_svg=true )   where T <: Number
    A, b = _forwardsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A[end:-1:1, end:-1:1], b[end:-1:1], var_name=var_name)
    lines = _relabel_cascade(lines, size(A, 1); var_name=var_name)
    if render_svg
        return _render_backsubst_svg(lines; fig_scale=fig_scale, tmp_dir=pb.tmp_dir, keep_file=pb.keep_file)
    end
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_forwardsubstitution!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x", fig_scale=1, render_svg=true )   where T <: Number"""
function show_forwardsubstitution!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x", fig_scale=1, render_svg=true )   where T <: Number
    A, b = _forwardsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A[end:-1:1, end:-1:1], b[end:-1:1], var_name=var_name)
    lines = _relabel_cascade(lines, size(A, 1); var_name=var_name)
    if render_svg
        return _render_backsubst_svg(lines; fig_scale=fig_scale, tmp_dir=pb.tmp_dir, keep_file=pb.keep_file)
    end
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_forwardsubstitution!(  pb::ShowGe{T}; b_col=1, var_name::String="x", fig_scale=1, render_svg=true )   where T <: Integer"""
function show_forwardsubstitution!(  pb::ShowGe{T}; b_col=1, var_name::String="x", fig_scale=1, render_svg=true )   where T <: Integer
    A, b = _forwardsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A[end:-1:1, end:-1:1], b[end:-1:1], var_name=var_name)
    lines = _relabel_cascade(lines, size(A, 1); var_name=var_name)
    if render_svg
        return _render_backsubst_svg(lines; fig_scale=fig_scale, tmp_dir=pb.tmp_dir, keep_file=pb.keep_file)
    end
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_solution!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x", fig\\_scale=1 )   where T <: Number"""
function show_solution!(  pb::ShowGe{Complex{Rational{T}}}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Number
    A, b = _backsub_ref(pb; b_col=b_col)
    tex = load_la_figures().standard_solution_tex(A, b, var_name=var_name)
    return _render_solution_svg(tex; fig_scale=fig_scale, tmp_dir=pb.tmp_dir)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_solution!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x", fig\\_scale=1 )   where T <: Number"""
function show_solution!(  pb::ShowGe{Rational{T}}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Number
    A, b = _backsub_ref(pb; b_col=b_col)
    tex = load_la_figures().standard_solution_tex(A, b, var_name=var_name)
    return _render_solution_svg(tex; fig_scale=fig_scale, tmp_dir=pb.tmp_dir)
end
# --------------------------------------------------------------------------------------------------------------
raw"""function show_solution!(  pb::ShowGe{T}; b_col=1, var_name::String="x", fig\\_scale=1 )   where T <: Integer"""
function show_solution!(  pb::ShowGe{T}; b_col=1, var_name::String="x", fig_scale=1 )   where T <: Integer
    A, b = _backsub_ref(pb; b_col=b_col)
    tex = load_la_figures().standard_solution_tex(A, b, var_name=var_name)
    return _render_solution_svg(tex; fig_scale=fig_scale, tmp_dir=pb.tmp_dir)
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
    show_forwardsubstitution(A, b; var_name="x", fig_scale=1, tmp_dir="tmp", keep_file=nothing, render_svg=true)

Render the forward-substitution cascade for the lower-triangular system `A * x = b`
using the la_figures backsubstitution cascade on a reversed system, then relabeling indices.
Supports Integer/Float as well as exact `Rational` and `Complex{Rational}` inputs
converted to tuples for exact SymPy reconstruction.
"""
function show_forwardsubstitution(A, b; var_name::String="x", fig_scale=1, tmp_dir="tmp", keep_file=nothing, render_svg=true)
    A2 = (A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(A) : A
    b2 = (b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(b) : b
    lines = load_la_figures().backsubstitution_tex(A2[end:-1:1, end:-1:1], b2[end:-1:1], var_name=var_name)
    lines = _relabel_cascade(lines, size(A, 1); var_name=var_name)
    if render_svg
        return _render_backsubst_svg(lines; fig_scale=fig_scale, tmp_dir=tmp_dir, keep_file=keep_file)
    end
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
    _ensure_pythoncall()
    return Base.invokelatest(PythonCall.pyconvert, String, s)
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

function _matrices_are_strings(mats)
    for row in mats
        for cell in row
            if cell === nothing
                continue
            end
            if cell isa AbstractArray
                for v in cell
                    if v === nothing
                        continue
                    end
                    return v isa AbstractString
                end
            else
                return cell isa AbstractString
            end
        end
    end
    return false
end

function _ge_block_to_list(block)
    if block === nothing || block === :none
        return nothing
    end
    if block isa AbstractArray
        rows = Vector{Any}()
        for i in axes(block, 1)
            row = Vector{Any}()
            for j in axes(block, 2)
                push!(row, block[i, j])
            end
            push!(rows, row)
        end
        return rows
    end
    return block
end

function _ge_grid_to_lists(mats)
    return [[_ge_block_to_list(block) for block in row] for row in mats]
end

function _ge_normalize_grid(mats)
    if mats isa AbstractVector
        if isempty(mats)
            return mats
        end
        first = mats[1]
        if first isa AbstractArray && !(first isa AbstractVector)
            return [[m] for m in mats]
        end
    end
    return mats
end

function _ge_to_pylist(obj)
    if obj isa AbstractArray
        _ensure_pythoncall()
        return Base.invokelatest(PythonCall.pylist, [_ge_to_pylist(x) for x in obj])
    end
    return obj
end

function _needs_shift_locs(obj)
    if obj isa Tuple && length(obj) == 2
        if all(x -> x isa Integer, obj)
            return any(x -> x == 0, obj)
        end
        if all(x -> x isa Tuple && length(x) == 2, obj)
            return _needs_shift_locs(obj[1]) || _needs_shift_locs(obj[2])
        end
    elseif obj isa AbstractVector
        return any(_needs_shift_locs, obj)
    end
    return false
end

function _shift_loc_pair(pair::Tuple)
    return (pair[1] + 1, pair[2] + 1)
end

function _shift_locs(obj)
    if obj isa Tuple && length(obj) == 2 && all(x -> x isa Integer, obj)
        return _shift_loc_pair(obj)
    elseif obj isa Tuple && length(obj) == 2 && all(x -> x isa Tuple && length(x) == 2, obj)
        return (_shift_locs(obj[1]), _shift_locs(obj[2]))
    elseif obj isa AbstractVector
        return [_shift_locs(x) for x in obj]
    end
    return obj
end

function _shift_specs(specs, loc_index::Int)
    if specs === nothing
        return nothing
    end
    if specs isa AbstractVector && !isempty(specs) && all(x -> x isa AbstractVector, specs)
        out = Vector{Any}(undef, length(specs))
        for (i, spec) in enumerate(specs)
            spec2 = copy(spec)
            if length(spec2) >= loc_index && _needs_shift_locs(spec2[loc_index])
                spec2[loc_index] = _shift_locs(spec2[loc_index])
            end
            out[i] = spec2
        end
        return out
    end
    if specs isa AbstractVector
        spec2 = copy(specs)
        if length(spec2) >= loc_index && _needs_shift_locs(spec2[loc_index])
            spec2[loc_index] = _shift_locs(spec2[loc_index])
        end
        return spec2
    end
    return specs
end

function matrixlayout_ge( matrices; Nrhs=0, formater=to_latex, pivot_list=nothing, bg_for_entries=nothing,
             variable_colors=["blue","black"], pivot_colors=["blue","yellow!40"], pivot_text_color=nothing,
             ref_path_list=nothing, comment_list=[], variable_summary=nothing, array_names=nothing,
             start_index=1, func=nothing, fig_scale=nothing, tmp_dir=nothing, keep_file=nothing, kwargs... )
    mats = matrices
    if !_matrices_are_strings(mats)
        mats = formater(mats)
    end
    mats = _ge_normalize_grid(mats)
    mats = _ge_grid_to_lists(mats)
    # Keep legacy 0-based coordinates for pivot/background specs; ge_convenience expects them.
    pivot_list = _ge_to_pylist(pivot_list)
    bg_for_entries = _ge_to_pylist(bg_for_entries)
    ref_path_list = _ge_to_pylist(ref_path_list)
    comment_list = _ge_to_pylist(comment_list)
    variable_summary = _ge_to_pylist(variable_summary)
    array_names = _ge_to_pylist(array_names)
    _ensure_pythoncall()
    builtins = _pyimport("builtins")
    py_str = Base.invokelatest(PythonCall.pygetattr, builtins, "str")
    ge_conv = _pyimport("la_figures.ge_convenience")
    ge_fn = Base.invokelatest(PythonCall.pygetattr, ge_conv, "ge")
    if pivot_text_color === nothing
        pivot_text_color = pivot_colors[1]
    end
    mats_py = _ge_to_pylist(mats)
    svg = _pycall(
        ge_fn,
        mats_py;
        Nrhs=Nrhs,
        formatter=py_str,
        pivot_list=pivot_list,
        bg_for_entries=bg_for_entries,
        variable_colors=variable_colors,
        pivot_text_color=pivot_text_color,
        ref_path_list=ref_path_list,
        comment_list=comment_list,
        variable_summary=variable_summary,
        array_names=array_names,
        start_index=start_index,
        func=func,
        fig_scale=fig_scale,
        tmp_dir=tmp_dir,
        keep_file=keep_file,
    )
    svg_str = Base.invokelatest(PythonCall.pyconvert, String, svg)
    return SVGOut(svg_str), nothing
end

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
