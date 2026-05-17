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

raw"""pb = ShowGE{T}(A::AbstractMatrix{T}; output_dir="/tmp/la/run", keep_file="/tmp/la/run/show\\_layout") where T <: Number
  <br>pb = ShowGE{T}(A::AbstractMatrix{T}; tmp_dir="/tmp/la/run", output_dir="/tmp/la/run", keep_file="/tmp/la/run/show\\_layout") where T <: Number
  <br>pb = ShowGE{T}(A::AbstractMatrix{T}, B::AbstractVector{T}; output_dir="/tmp/la/run", keep_file="/tmp/la/run/show\\_layout") where T <: Number
  <br>pb = ShowGE{T}(A::AbstractMatrix{T}, B::AbstractMatrix{T}; output_dir="/tmp/la/run", keep_file="/tmp/la/run/show\\_layout") where T <: Number
  <br>pb = ShowGE{Rational{T}}(A::AbstractMatrix{T}, B::AbstractMatrix{T}; output_dir="/tmp/la/run", keep_file="/tmp/la/run/show\\_layout") where T <: Number
  <br>pb = ShowGE{Rational{T}}(A::AbstractMatrix{T}, B::AbstractVector{T}; output_dir="/tmp/la/run", keep_file="/tmp/la/run/show\\_layout") where T <: Number
  <br>pb = ShowGE{Complex{Rational{T}}}(A::AbstractMatrix{Complex{T}}; output_dir="/tmp/la/run", keep_file="/tmp/la/run/show\\_layout") where T <: Number
  <br>pb = ShowGE{Rational{T}}(A::AbstractMatrix{T}; output_dir="/tmp/la/run", keep_file="/tmp/la/run/show\\_layout") where T <: Number
  <br>pb = ShowGE{Complex{Rational{T}}}(A::AbstractMatrix{Complex{T}}, B::AbstractVector{Complex{T}}; output_dir="/tmp/la/run", keep_file="/tmp/la/run/show\\_layout") where T <: Number
  <br>pb = ShowGE{Complex{Rational{T}}}(A::AbstractMatrix{Complex{T}}, B::AbstractMatrix{Complex{T}}; output_dir="/tmp/la/run", keep_file="/tmp/la/run/show\\_layout") where T <: Number"""
mutable struct ShowGE{T<:Number}
    tmp_dir
    keep_file
    A
    B
    n_rhs
    normal_eq::Bool

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
    xp
    xh
    status
    rhs_status
    rhs_consistent


  function ShowGE(A::AbstractMatrix; tmp_dir="/tmp/la/run", output_dir=nothing, keep_file="/tmp/la/run/show_layout")
	  ShowGE{eltype(A)}(A; tmp_dir=tmp_dir, output_dir=output_dir, keep_file=keep_file)
  end
  function ShowGE(A::AbstractMatrix, b; tmp_dir="/tmp/la/run", output_dir=nothing, keep_file="/tmp/la/run/show_layout")
    ShowGE{eltype(A)}(A, b; tmp_dir=tmp_dir, output_dir=output_dir, keep_file=keep_file)
  end
  function ShowGE{T}(A::AbstractMatrix{T}; tmp_dir="/tmp/la/run", output_dir=nothing, keep_file="/tmp/la/run/show_layout") where T <: Number
      new(_resolve_output_dir(output_dir, tmp_dir), keep_file, A, nothing, 0, false)
  end
  function ShowGE{T}(A::AbstractMatrix{T}, B::AbstractVector{T}; tmp_dir="/tmp/la/run", output_dir=nothing, keep_file="/tmp/la/run/show_layout") where T <: Number
      new(_resolve_output_dir(output_dir, tmp_dir), keep_file, A, B, size(B,2), false)
  end
  function ShowGE{T}(A::AbstractMatrix{T}, B::AbstractMatrix{T}; tmp_dir="/tmp/la/run", output_dir=nothing, keep_file="/tmp/la/run/show_layout") where T <: Number
      new(_resolve_output_dir(output_dir, tmp_dir), keep_file, A, B, size(B,2), false)
  end

  function ShowGE{Rational{T}}(A::AbstractMatrix{T}; tmp_dir="/tmp/la/run", output_dir=nothing, keep_file="/tmp/la/run/show_layout") where T <: Number
      new(_resolve_output_dir(output_dir, tmp_dir), keep_file, Rational{T}.(A), nothing, 0, false)
  end
  function ShowGE{Rational{T}}(A::AbstractMatrix{T}, B::AbstractVector{T}; tmp_dir="/tmp/la/run", output_dir=nothing, keep_file="/tmp/la/run/show_layout") where T <: Number
      new(_resolve_output_dir(output_dir, tmp_dir), keep_file, Rational{T}.(A),Rational{T}.(B),size(B,2), false)
  end
  function ShowGE{Rational{T}}(A::AbstractMatrix{T}, B::AbstractMatrix{T}; tmp_dir="/tmp/la/run", output_dir=nothing, keep_file="/tmp/la/run/show_layout") where T <: Number
      new(_resolve_output_dir(output_dir, tmp_dir), keep_file, Rational{T}.(A),Rational{T}.(B),size(B,2), false)
  end

  function ShowGE{Complex{Rational{T}}}(A::AbstractMatrix{Complex{T}}; tmp_dir="/tmp/la/run", output_dir=nothing, keep_file="/tmp/la/run/show_layout") where T <: Number
    new(_resolve_output_dir(output_dir, tmp_dir), keep_file, Complex{Rational{T}}.(A), nothing, 0, false)
  end
  function ShowGE{Complex{Rational{T}}}(A::AbstractMatrix{Complex{T}}, B::AbstractVector{Complex{T}}; tmp_dir="/tmp/la/run", output_dir=nothing, keep_file="/tmp/la/run/show_layout") where T <: Number
    new(_resolve_output_dir(output_dir, tmp_dir), keep_file, Complex{Rational{T}}.(A),Complex{Rational{T}}.(B),size(B,2), false)
  end
  function ShowGE{Complex{Rational{T}}}(A::AbstractMatrix{Complex{T}}, B::AbstractMatrix{Complex{T}}; tmp_dir="/tmp/la/run", output_dir=nothing, keep_file="/tmp/la/run/show_layout") where T <: Number
      new(_resolve_output_dir(output_dir, tmp_dir), keep_file, Complex{Rational{T}}.(A),Complex{Rational{T}}.(B),size(B,2), false)
  end
end

# --------------------------------------------------------------------------------------------------------------
"""
    ref!(pb; n_rhs=:None, gj=false, normal_eq=false)

Compute REF/RREF data for a `ShowGE` problem and attach pivot metadata.
"""
function ref!( pb::ShowGE{T}; n_rhs=:None, gj::Bool=false, normal_eq::Bool=false )  where T <: Number
    M,N = size(pb.A)
    if pb.B !== nothing
        A = [pb.A pb.B]
        if n_rhs != :None
            pb.n_rhs = n_rhs
        end
    else
       A = pb.A
       pb.n_rhs = 0
    end
    if normal_eq
      pb.matrices, pb.pivot_cols, pb.desc = normal_eq_reduce_to_ref( A, n=N, gj=gj );
      sz = (N,N)
    else
      pb.matrices, pb.pivot_cols, pb.desc = reduce_to_ref( A, n=N, gj=gj );
      sz = (M,N)
    end
    pb.normal_eq = normal_eq
    pb.free_cols = filter(x -> !(x in pb.pivot_cols), 1:N)

    pb.pivot_list, pb.bg_for_entries, pb.ref_path_list, pb.basic_var = decorate_ge(pb.desc,pb.pivot_cols,sz; pivot_color="yellow!40");
    pb.rank = length( pb.pivot_cols )
    _compute_rhs_status!(pb)
    nothing
end
# --------------------------------------------------------------------------------------------------------------
function _rhs_val_to_tex(val)
    if val isa Rational
        return _rational_str(val)
    elseif val isa Complex{<:Rational}
        return _complex_rational_str(val)
    end
    return string(val)
end

function _format_rhs_label(rhs_labels::Vector{String})
    if isempty(rhs_labels)
        return ""
    end
    if length(rhs_labels) > 1
        return "\\left( " * join(rhs_labels, " \\mid ") * " \\right)"
    end
    return rhs_labels[1]
end

function _normal_eq_name_specs(n_rows::Int, rhs_labels::Vector{String})
    specs = Vector{Tuple{Tuple{Int,Int},String,String}}()
    if n_rows <= 0
        return specs
    end
    rhs0 = _format_rhs_label(rhs_labels)
    push!(specs, ((0, 1), "ar", "\$\\mathbf{ $(rhs0) }\$"))
    if n_rows <= 1
        return specs
    end
    rhs1_labels = ["A^T " * lbl for lbl in rhs_labels]
    rhs1 = _format_rhs_label(rhs1_labels)
    push!(specs, ((1, 0), "al", "\$\\mathbf{ A^T }\$"))
    push!(specs, ((1, 1), "ar", "\$\\mathbf{ $(rhs1) }\$"))
    for i in 2:(n_rows - 1)
        prod = join(["E_$(k)" for k in (i - 1):-1:1], " ")
        rhs_i_labels = [prod == "" ? lbl : "$(prod) $(lbl)" for lbl in rhs1_labels]
        rhs_i = _format_rhs_label(rhs_i_labels)
        push!(specs, ((i, 0), "al", "\$\\mathbf{ E_$(i - 1) }\$"))
        push!(specs, ((i, 1), "ar", "\$\\mathbf{ $(rhs_i) }\$"))
    end
    return specs
end

function _inconsistent_rhs_value(pb::ShowGE{T}, b_col::Int) where T <: Number
    Ab = pb.matrices[end][end]
    n = size(pb.A, 2)
    m = size(Ab, 1)
    for i in 1:m
        row_zero = true
        for j in 1:n
            if Ab[i, j] != 0
                row_zero = false
                break
            end
        end
        if row_zero && Ab[i, n + b_col] != 0
            return Ab[i, n + b_col]
        end
    end
    return nothing
end

function _rhs_col_count(pb::ShowGE)
    if pb.B === nothing
        return 0
    end
    return size(pb.B, 2)
end

function _compute_rhs_status!(pb::ShowGE{T}) where T <: Number
    n_rhs = _rhs_col_count(pb)
    if n_rhs == 0
        pb.status = :none
        pb.rhs_status = Symbol[]
        pb.rhs_consistent = Bool[]
        return pb
    end
    Ab = pb.matrices[end][end]
    n = size(pb.A, 2)
    m = size(Ab, 1)
    rhs_status = Vector{Symbol}(undef, n_rhs)
    rhs_consistent = Vector{Bool}(undef, n_rhs)
    for k in 1:n_rhs
        status = :consistent
        for i in 1:m
            row_zero = true
            for j in 1:n
                if Ab[i, j] != 0
                    row_zero = false
                    break
                end
            end
            if row_zero && Ab[i, n + k] != 0
                status = :inconsistent
                break
            end
        end
        rhs_status[k] = status
        rhs_consistent[k] = (status == :consistent)
    end
    pb.rhs_status = rhs_status
    pb.rhs_consistent = rhs_consistent
    if all(s -> s == :consistent, rhs_status)
        pb.status = :consistent
    elseif all(s -> s == :inconsistent, rhs_status)
        pb.status = :inconsistent
    else
        pb.status = :mixed
    end
    return pb
end
# --------------------------------------------------------------------------------------------------------------
"""
    show_layout!(pb; array_names=nothing, show_variables=true, fig_scale=1, output_dir=nothing, output_stem=nothing)

Render the GE layout for a `ShowGE` problem as SVG.
Use `output_dir` and `output_stem` to retain TeX/PDF/SVG artifacts.
"""
function _keep_file_output_parts(keep_file)
    if keep_file === nothing
        return nothing, nothing
    end
    path = String(keep_file)
    dir = dirname(path)
    stem = splitext(basename(path))[1]
    return dir, stem
end

function show_layout!(  pb::ShowGE{T}; array_names=nothing, show_variables=true, fig_scale=1, output_dir=nothing, output_stem=nothing, render_opts=nothing )   where T <: Number
    n_rhs = _rhs_col_count(pb)
    if array_names === nothing
        if n_rhs == 0
            array_names = ["E", "A"]
        elseif n_rhs == 1
            array_names = ["E", ["A", "b"]]
        else
            array_names = ["E", ["A", "B"]]
        end
    end
    if pb.normal_eq
        if !(array_names isa AbstractDict)
            rhs_labels = String[]
            try
                _, rhs = array_names
                rhs_labels = [string(x) for x in rhs]
            catch
                rhs_labels = ["A"]
            end
            name_specs = _normal_eq_name_specs(length(pb.matrices), rhs_labels)
            array_names = Dict("name_specs" => name_specs)
        end
            rhs_status = isdefined(pb, :rhs_status) ? [string(s) for s in pb.rhs_status] : nothing
            keep_dir, keep_stem = _keep_file_output_parts(pb.keep_file)
            resolved_output_dir = output_dir !== nothing ? output_dir : (keep_dir !== nothing ? keep_dir : pb.tmp_dir)
            resolved_output_stem = output_stem !== nothing ? output_stem : keep_stem
            svg = matrixlayout_ge(
                pb.matrices;
                n_rhs=pb.n_rhs,
            pivot_list=pb.pivot_list,
            bg_for_entries=pb.bg_for_entries,
            ref_path_list=pb.ref_path_list,
            variable_summary=show_variables ? pb.basic_var : nothing,
            variable_colors=["red", "black"],
                rhs_status=rhs_status,
                array_names=array_names,
                fig_scale=fig_scale,
                output_dir=resolved_output_dir,
                output_stem=resolved_output_stem,
                render_opts=render_opts,
            )
        pb.h = svg
        return svg
    end
    if isdefined(pb, :matrices) && pb.matrices !== nothing && length(pb.matrices) > 1
            rhs_status = isdefined(pb, :rhs_status) ? [string(s) for s in pb.rhs_status] : nothing
            keep_dir, keep_stem = _keep_file_output_parts(pb.keep_file)
            resolved_output_dir = output_dir !== nothing ? output_dir : (keep_dir !== nothing ? keep_dir : pb.tmp_dir)
            resolved_output_stem = output_stem !== nothing ? output_stem : keep_stem
            svg = matrixlayout_ge(
                pb.matrices;
                n_rhs=pb.n_rhs,
            pivot_list=pb.pivot_list,
            bg_for_entries=pb.bg_for_entries,
            ref_path_list=pb.ref_path_list,
            variable_summary=show_variables ? pb.basic_var : nothing,
            variable_colors=["red", "black"],
                rhs_status=rhs_status,
                array_names=array_names,
                fig_scale=fig_scale,
                output_dir=resolved_output_dir,
                output_stem=resolved_output_stem,
                render_opts=render_opts,
            )
        pb.h = svg
        return svg
    end
    la = load_la_figures()
    rhs = isdefined(pb, :B) ? pb.B : nothing
    ge_tbl_svg = _pygetattr_fallback(la, :ge_tbl_svg, "la_figures.ge_convenience")
    rhs_status = isdefined(pb, :rhs_status) ? pb.rhs_status : Symbol[]
    rhs_status_str = [string(s) for s in rhs_status]
    keep_dir, keep_stem = _keep_file_output_parts(pb.keep_file)
    resolved_output_dir = output_dir !== nothing ? output_dir : (keep_dir !== nothing ? keep_dir : pb.tmp_dir)
    pb.h = _pycall(ge_tbl_svg, pb.A, rhs;
        show_pivots=true,
        fig_scale=fig_scale,
        variable_summary=show_variables ? pb.basic_var : nothing,
        variable_colors=["red", "black"],
        rhs_status=rhs_status_str,
        array_names=array_names,
        output_dir=resolved_output_dir,
        render_opts=render_opts,
    )
    _ensure_pythoncall()
    svg_str = Base.invokelatest(PythonCall.pyconvert, String, pb.h)
    svg = SVGOut(svg_str)
    pb.h = svg
    return svg
end

# --------------------------------------------------------------------------------------------------------------
"""
    show_system(pb; b_col=1, var_name="x", fig_scale=1)

Render the linear system associated with a `ShowGE` problem.
"""
function _system_matrix_rhs(pb::ShowGE{T}; b_col=1) where T <: Number
    A = pb.A
    if isdefined(pb, :B) && b_col isa Integer && 1 <= b_col <= size(pb.B, 2)
       b = pb.B[:,b_col]
    else
       b = zeros(eltype(pb.A), size(A,1), 1)
    end
    if A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}
        A = _encode_exact.(A)
    end
    if b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}
        b = _encode_exact.(b)
    end
    return A, b
end

_resolve_output_dir(output_dir, tmp_dir, fallback=nothing) = output_dir !== nothing ? output_dir : (tmp_dir !== nothing ? tmp_dir : fallback)

function show_system(  pb::ShowGE{T}; b_col=1, var_name::String="x", fig_scale=1, output_dir=nothing, tmp_dir=nothing, render_opts=nothing) where T <: Number
    A, b = _system_matrix_rhs(pb; b_col=b_col)
    la = load_la_figures()
    linear_system_tex = _pygetattr(la, :linear_system_tex)
    tex = _pycall(linear_system_tex, A, b; var_name=var_name)
    _ensure_pythoncall()
    tex = Base.invokelatest(PythonCall.pyconvert, String, tex)
    bs = _pyimport("matrixlayout.backsubst")
    backsubst_svg = _pygetattr(bs, :backsubst_svg)
    resolved_output_dir = _resolve_output_dir(output_dir, tmp_dir, pb.tmp_dir)
    svg = _pycall(backsubst_svg; system_txt=tex, show_system=true,
                  show_cascade=false, show_solution=false,
                  fig_scale=fig_scale, output_dir=resolved_output_dir,
                  render_opts=render_opts)
    return _show_svg(svg)
end
"""
    create_cascade!(pb::ShowGE; b_col=1, var_name="x")

Initialize cascade state for substitution displays.
"""
function create_cascade!(  pb::ShowGE{T}; b_col=1, var_name::String="x" ) where T <: Number
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

function _backsub_ref(pb::ShowGE; b_col=1)
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

function _forwardsub_ref(pb::ShowGE; b_col=1)
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

function _render_backsubst_svg(lines; fig_scale=nothing, output_dir=nothing, render_opts=nothing)
    ml = load_matrixlayout()
    backsubst_svg = _pygetattr(ml, :backsubst_svg)
    kwargs = Dict{Symbol, Any}()
    render_opts = _normalize_render_opts(render_opts; tmp_dir=output_dir)
    kwargs[:cascade_txt] = _ge_to_pylist(lines)
    kwargs[:show_system] = false
    kwargs[:show_cascade] = true
    kwargs[:show_solution] = false
    if fig_scale !== nothing
        kwargs[:fig_scale] = fig_scale
    end
    if output_dir !== nothing
        kwargs[:output_dir] = output_dir
    end
    svg = _pycall(backsubst_svg; kwargs..., render_opts=render_opts)
    return _show_svg(svg)
end

function _render_solution_svg(solution_tex; fig_scale=nothing, output_dir=nothing, render_opts=nothing)
    ml = load_matrixlayout()
    backsubst_svg = _pygetattr(ml, :backsubst_svg)
    kwargs = Dict{Symbol, Any}()
    render_opts = _normalize_render_opts(render_opts; tmp_dir=output_dir)
    kwargs[:solution_txt] = solution_tex
    kwargs[:show_system] = false
    kwargs[:show_cascade] = false
    kwargs[:show_solution] = true
    if fig_scale !== nothing
        kwargs[:fig_scale] = fig_scale
    end
    if output_dir !== nothing
        kwargs[:output_dir] = output_dir
    end
    svg = _pycall(backsubst_svg; kwargs..., render_opts=render_opts)
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
"""
    show_backsubstitution!(pb; b_col=1, var_name="x", fig_scale=1)

Render the back-substitution cascade for a `ShowGE` problem.
"""
function show_backsubstitution!(  pb::ShowGE{T}; b_col=1, var_name::String="x", fig_scale=1, render_opts=nothing ) where T <: Number
    if isdefined(pb, :rhs_status) && b_col isa Integer && b_col <= length(pb.rhs_status) && pb.rhs_status[b_col] == :inconsistent
        val = _inconsistent_rhs_value(pb, b_col)
        rhs_txt = val === nothing ? "?" : _rhs_val_to_tex(val)
        lines = [string("0 = ", rhs_txt), "\\text{No Solution}"]
        return _render_backsubst_svg(lines; fig_scale=fig_scale, output_dir=pb.tmp_dir, render_opts=render_opts)
    end
    A, b = _backsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A, b, var_name=var_name)
    return _render_backsubst_svg(lines; fig_scale=fig_scale, output_dir=pb.tmp_dir, render_opts=render_opts)
end
# --------------------------------------------------------------------------------------------------------------
function show_forwardsubstitution!(  pb::ShowGE{T}; b_col=1, var_name::String="x", fig_scale=1, render_svg=true, render_opts=nothing ) where T <: Number
    A, b = _forwardsub_ref(pb; b_col=b_col)
    lines = load_la_figures().backsubstitution_tex(A[end:-1:1, end:-1:1], b[end:-1:1], var_name=var_name)
    lines = _relabel_cascade(lines, size(A, 1); var_name=var_name)
    if render_svg
        return _render_backsubst_svg(lines; fig_scale=fig_scale, output_dir=pb.tmp_dir, render_opts=render_opts)
    end
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
"""
    show_solution!(pb; b_col=1, var_name="x", fig_scale=1)

Render the solution vector/form for a `ShowGE` problem.
"""
function show_solution!(  pb::ShowGE{T}; b_col=1, var_name::String="x", fig_scale=1, render_opts=nothing ) where T <: Number
    if isdefined(pb, :rhs_status) && b_col isa Integer && b_col <= length(pb.rhs_status) && pb.rhs_status[b_col] == :inconsistent
        return Vector{T}()
    end
    A, b = _backsub_ref(pb; b_col=b_col)
    tex = load_la_figures().standard_solution_tex(A, b, var_name=var_name)
    return _render_solution_svg(tex; fig_scale=fig_scale, output_dir=pb.tmp_dir, render_opts=render_opts)
end
# ==============================================================================================================
raw"""
    show_backsubstitution(A, b; var_name="x", fig_scale=1, output_dir=nothing, tmp_dir="/tmp/la/run", render_svg=true)

    Render the back-substitution cascade for the upper-triangular system `A * x = b`
    using `la_figures.backsubstitution_tex`. Works with Integer/Float as well as
    exact `Rational` and `Complex{Rational}` inputs (those are converted to tuples so
    SymPy reconstructs exact rationals on the Python side).
"""
function show_backsubstitution(A, b; var_name::String="x", fig_scale=1, output_dir=nothing, tmp_dir="/tmp/la/run", render_svg=true, render_opts=nothing)
    A2 = (A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(A) : A
    b2 = (b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(b) : b
    lines = load_la_figures().backsubstitution_tex(A2, b2, var_name=var_name)
    if render_svg
        return _render_backsubst_svg(lines; fig_scale=fig_scale, output_dir=_resolve_output_dir(output_dir, tmp_dir), render_opts=render_opts)
    end
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
raw"""
    show_forwardsubstitution(A, b; var_name="x", fig_scale=1, output_dir=nothing, tmp_dir="/tmp/la/run", render_svg=true)

Render the forward-substitution cascade for the lower-triangular system `A * x = b`
using the la_figures backsubstitution cascade on a reversed system, then relabeling indices.
Supports Integer/Float as well as exact `Rational` and `Complex{Rational}` inputs
converted to tuples for exact SymPy reconstruction.
"""
function show_forwardsubstitution(A, b; var_name::String="x", fig_scale=1, output_dir=nothing, tmp_dir="/tmp/la/run", render_svg=true, render_opts=nothing)
    A2 = (A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(A) : A
    b2 = (b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(b) : b
    lines = load_la_figures().backsubstitution_tex(A2[end:-1:1, end:-1:1], b2[end:-1:1], var_name=var_name)
    lines = _relabel_cascade(lines, size(A, 1); var_name=var_name)
    if render_svg
        return _render_backsubst_svg(lines; fig_scale=fig_scale, output_dir=_resolve_output_dir(output_dir, tmp_dir), render_opts=render_opts)
    end
    return _display_cascade(lines)
end
# ==============================================================================================================
"""
    solutions(pb::ShowGE) -> xp, xh
"""
function solutions(pb::ShowGE{T} )   where T <: Number
    _, N = size(pb.A)
    matrices, pivot_cols, _ = reduce_to_ref(pb.matrices[end][end][1:pb.rank, 1:end], n=N, gj=true)

    if _rhs_col_count(pb) > 0
        rhs_consistent = isdefined(pb, :rhs_consistent) ? pb.rhs_consistent : Bool[]
        keep = [i for (i, ok) in enumerate(rhs_consistent) if ok]
        if isempty(keep)
            Xp = zeros(T, N, 0)
        else
            Xp = zeros(T, N, length(keep))
            F = matrices[end][end][1:pb.rank, N+1:end]
            Xp[pivot_cols, :] = F[:, keep]
        end
    else
        Xp = zeros(T, N, 1)
    end

    if length(pb.free_cols) > 0
        Xh = zeros(T, N, N - pb.rank)
        F = matrices[end][end][1:pb.rank, pb.free_cols]
        for (col, row) in enumerate(pb.free_cols)
            Xh[row, col] = 1
        end
        Xh[pivot_cols, :] = -F
    else
        Xh = zeros(T, N, 1)
    end
    if size(Xh, 2) == 0 || all(x -> x == 0, Xh)
        Xh = zeros(T, N, 1)
    end
    pb.xp = Xp
    pb.xh = Xh
    return Xp, Xh
end
# ------------------------------------------------------------------------------------------
function solve!(pb::ShowGE{T} )   where T <: Number
    pb.xp, pb.xh = solutions(pb)
end
# ==============================================================================================================
"""
    rhs_block(pb::ShowGE; step=:final, b_col=nothing)

Return the RHS block from a GE stack. Use `step` to select a specific matrix
stack index (1-based) or `:final` for the last step. If `b_col` is provided,
returns a single RHS column vector; otherwise returns the full RHS block.
"""
function rhs_block(pb::ShowGE{T}; step=:final, b_col=nothing) where T <: Number
    if _rhs_col_count(pb) == 0
        return nothing
    end
    mats = pb.matrices
    idx = step === :final ? length(mats) : step
    Ab = mats[idx][end]
    n = size(pb.A, 2)
    rhs = Ab[:, n + 1:end]
    if b_col === nothing
        return rhs
    end
    return rhs[:, b_col]
end
# ==============================================================================================================
# function column_view( Xp, Xh, pivot_cols, rhs )
# end
# ==============================================================================================================
#function homogeneous_solution(pb::ShowGE{Complex{Rational{T}}}; b_col=1 )   where T <: Number)
#  N = size(pb.A,2)
#  matrices, pivot_cols, desc = reduce_to_ref( pb.matrices[end][end][:,1:N], n=N, gj=true );
#  Xh = similar(pb.A, size(pb.A,1), A - pb.rank)
#end
# ==============================================================================================================
raw"""function show_ge_final( matrices, desc, pivot_cols; n_rhs=0, formatter=to_latex, pivot_list=nothing, bg_for_entries=nothing, <br>
             variable_colors=["blue","black"], pivot_colors=["blue","yellow!40"],  <br>
             ref_path_list=nothing, comment_list=[], variable_summary=nothing, array_names=nothing, <br>
             start_index=1, func=nothing, fig_scale=nothing, output_dir=nothing, output_stem=nothing, tmp_dir=nothing )
 """
function julia_ge( matrices, desc, pivot_cols; n_rhs=0, formatter=to_latex, pivot_list=nothing, bg_for_entries=nothing,
             variable_colors=["blue","black"], pivot_colors=["blue","yellow!40"],
             ref_path_list=nothing, comment_list=[], variable_summary=nothing, array_names=nothing,
             start_index=1, func=nothing, fig_scale=nothing, output_dir=nothing, output_stem=nothing, tmp_dir=nothing,
             render_opts=nothing )
    Ab = matrices[end][end]
    nrhs = n_rhs isa AbstractArray ? sum(n_rhs) : n_rhs
    if nrhs > 0
        A = Ab[:, 1:(end - nrhs)]
        rhs = Ab[:, (end - nrhs + 1):end]
    else
        A = Ab
        rhs = nothing
    end
    la = load_la_figures()
    ge_tbl_svg = _pygetattr_fallback(la, :ge_tbl_svg, "la_figures.ge_convenience")
    local_render_opts = render_opts === nothing ? Dict{String, Any}() : Dict{String, Any}(render_opts)
    resolved_output_dir = output_dir !== nothing ? output_dir : tmp_dir
    if !haskey(local_render_opts, "output_dir") && !haskey(local_render_opts, :output_dir)
        if output_dir !== nothing
            local_render_opts["output_dir"] = resolved_output_dir
        elseif tmp_dir !== nothing
            local_render_opts["output_dir"] = mktempdir(resolved_output_dir)
        end
    end
    call_kwargs = Dict{Symbol, Any}(
        :fig_scale => fig_scale,
        :array_names => array_names,
        :variable_summary => variable_summary,
        :variable_colors => variable_colors,
        :render_opts => local_render_opts,
    )
    if output_dir !== nothing
        call_kwargs[:output_dir] = resolved_output_dir
    end
    if output_stem !== nothing
        call_kwargs[:output_stem] = output_stem
    end
    s = _pycall(ge_tbl_svg, A, rhs; call_kwargs...)
    _ensure_pythoncall()
    return Base.invokelatest(PythonCall.pyconvert, String, s)
end
"""
    SVGOut(svg::String)

Wrapper type for SVG output from `show_ge_final`, enabling rich display in IJulia.
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
    show_ge_final(args...; kwargs...) -> SVGOut

Render only the final GE table from `matrices[end][end]`.
"""
show_ge_final(args...; kwargs...) = SVGOut(julia_ge(args...; kwargs...))

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
    if block isa LinearAlgebra.Adjoint
        parent_block = parent(block)
        # For non-numeric entries (e.g., latex string cells), avoid element-wise
        # adjoint application and perform a plain transpose of the container.
        block = eltype(parent_block) <: Number ? Matrix(block) : permutedims(parent_block)
    elseif block isa LinearAlgebra.Transpose
        block = permutedims(parent(block))
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
    if obj isa AbstractDict
        _ensure_pythoncall()
        return Base.invokelatest(PythonCall.pydict, Dict(k => _ge_to_pylist(v) for (k, v) in obj))
    end
    if obj isa AbstractArray
        _ensure_pythoncall()
        return Base.invokelatest(PythonCall.pylist, [_ge_to_pylist(x) for x in obj])
    end
    return obj
end

function _normalize_bg_specs(specs)
    if !(specs isa AbstractVector)
        return [specs]
    end
    if isempty(specs)
        return specs
    end
    if length(specs) >= 3 && !(specs[1] isa AbstractVector) && !(specs[2] isa AbstractVector)
        return [specs]
    end
    if !all(x -> x isa AbstractVector, specs)
        return [specs]
    end
    return specs
end

function _matrix_shape_from_grid(mats_raw, gM::Int, gN::Int)
    if mats_raw === nothing
        return 0, 0
    end
    try
        grid = mats_raw[gM + 1]
        mat = grid[gN + 1]
        if mat === nothing || mat isa AbstractString
            return 0, 0
        elseif mat isa AbstractMatrix
            return size(mat)
        elseif mat isa AbstractArray && ndims(mat) == 2
            return size(mat)
        elseif mat isa AbstractVector
            nrows = length(mat)
            ncols = nrows == 0 ? 0 : length(mat[1])
            return nrows, ncols
        end
    catch
    end
    return 0, 0
end

function _bg_for_entries_to_decorators(bg_for_entries, mats_raw=nothing)
    if bg_for_entries === nothing
        return nothing
    end
    _ensure_pythoncall()
    fmt = _pyimport("matrixlayout.formatting")
    make_decorator = Base.invokelatest(PythonCall.pygetattr, fmt, "make_decorator")
    sel_entry = Base.invokelatest(PythonCall.pygetattr, fmt, "sel_entry")
    sel_box = Base.invokelatest(PythonCall.pygetattr, fmt, "sel_box")
    specs = _normalize_bg_specs(bg_for_entries)
    decorators = Vector{Any}()
    for spec in specs
        if !(spec isa AbstractVector) || length(spec) < 3
            continue
        end
        gM = Int(spec[1])
        gN = Int(spec[2])
        nrows, ncols = _matrix_shape_from_grid(mats_raw, gM, gN)
        entries = spec[3]
        color = length(spec) >= 4 ? string(spec[4]) : "red!15"
        if !(entries isa AbstractVector)
            entries = [entries]
        end
        entry_selectors = Vector{Any}()
        decorator = Base.invokelatest(make_decorator; bg_color=color)
        for entry in entries
            if (entry isa AbstractVector || entry isa Tuple) && length(entry) == 2 &&
               (entry[1] isa AbstractVector || entry[1] isa Tuple) &&
               (entry[2] isa AbstractVector || entry[2] isa Tuple)
                i0, j0 = entry[1]
                i1, j1 = entry[2]
                push!(entry_selectors, Base.invokelatest(sel_box, (Int(i0), Int(j0)), (Int(i1), Int(j1))))
            else
                i0, j0 = entry
                push!(entry_selectors, Base.invokelatest(sel_entry, Int(i0), Int(j0)))
            end
        end
        push!(decorators, Dict(
            "grid" => (gM, gN),
            "entries" => entry_selectors,
            "decorator" => decorator,
        ))
    end
    return isempty(decorators) ? nothing : decorators
end

function _prepare_ge_mats(matrices, formatter)
    mats = matrices
    if !_matrices_are_strings(mats)
        mats = formatter(mats)
    end
    mats = _ge_normalize_grid(mats)
    return _ge_grid_to_lists(mats)
end

function _merge_bg_decorators(bg_for_entries, decorators, mats)
    decorators_from_bg = _bg_for_entries_to_decorators(bg_for_entries, mats)
    if decorators !== nothing && !(decorators isa AbstractVector)
        decorators = [decorators]
    end
    if decorators_from_bg !== nothing
        bg_for_entries = nothing
        decorators = decorators === nothing ? decorators_from_bg : vcat(decorators_from_bg, decorators)
    end
    return bg_for_entries, decorators
end

function _call_ge_convenience(mats_py; kwargs...)
    _ensure_pythoncall()
    ge_conv = _pyimport("la_figures.ge_convenience")
    ge_fn = Base.invokelatest(PythonCall.pygetattr, ge_conv, "ge")
    return _pycall(ge_fn, mats_py; kwargs...)
end

function _ge_pyify_payload(mats, pivot_list, bg_for_entries, ref_path_list, comment_list, variable_summary, rhs_status, array_names, specs, decorators)
    return (
        mats=_ge_to_pylist(mats),
        pivot_list=_ge_to_pylist(pivot_list),
        bg_for_entries=_ge_to_pylist(bg_for_entries),
        ref_path_list=_ge_to_pylist(ref_path_list),
        comment_list=_ge_to_pylist(comment_list),
        variable_summary=_ge_to_pylist(variable_summary),
        rhs_status=_ge_to_pylist(rhs_status),
        array_names=_ge_to_pylist(array_names),
        specs=_ge_to_pylist(specs),
        decorators=_ge_to_pylist(decorators),
    )
end

function matrixlayout_ge( matrices; n_rhs=0, formatter=to_latex, pivot_list=nothing, bg_for_entries=nothing,
             variable_colors=["blue","black"], pivot_colors=["blue","yellow!40"], pivot_text_color=nothing,
             ref_path_list=nothing, comment_list=[], variable_summary=nothing, array_names=nothing,
             start_index=1, func=nothing, fig_scale=nothing, output_dir=nothing, output_stem=nothing, tmp_dir=nothing,
             render_opts=nothing, rhs_status=nothing, kwargs... )
    mats = _prepare_ge_mats(matrices, formatter)
    # Use 0-based coordinates for pivot/background specs; ge_convenience expects them.
    decorators = get(kwargs, :decorators, nothing)
    had_bg = bg_for_entries !== nothing
    bg_for_entries, decorators = _merge_bg_decorators(bg_for_entries, decorators, mats)
    # create_medium_nodes/create_extra_nodes are handled by la_figures.ge_convenience
    specs = get(kwargs, :specs, nothing)
    payload = _ge_pyify_payload(
        mats, pivot_list, bg_for_entries, ref_path_list,
        comment_list, variable_summary, rhs_status, array_names, specs, decorators
    )

    _ensure_pythoncall()
    builtins = _pyimport("builtins")
    py_str = Base.invokelatest(PythonCall.pygetattr, builtins, "str")
    nrhs_arg = n_rhs isa AbstractVector ? _ge_to_pylist(n_rhs) : n_rhs
    if pivot_text_color === nothing
        pivot_text_color = pivot_colors[1]
    end
    resolved_output_dir = output_dir !== nothing ? output_dir : tmp_dir
    svg = _call_ge_convenience(
        payload.mats;
        n_rhs=nrhs_arg,
        formatter=py_str,
        pivot_list=payload.pivot_list,
        bg_for_entries=payload.bg_for_entries,
        variable_colors=variable_colors,
        pivot_text_color=pivot_text_color,
        ref_path_list=payload.ref_path_list,
        comment_list=payload.comment_list,
        variable_summary=payload.variable_summary,
        rhs_status=payload.rhs_status,
        array_names=payload.array_names,
        specs=payload.specs,
        decorators=payload.decorators,
        start_index=start_index,
        func=func,
        fig_scale=fig_scale,
        output_dir=resolved_output_dir,
        output_stem=output_stem,
        render_opts=render_opts,
    )
    svg_str = Base.invokelatest(PythonCall.pyconvert, String, svg)
    return SVGOut(svg_str)
end

# ------------------------------------------------------------------------------------------
"""
    show_solution(matrices; var_name="x", fig_scale=1, output_dir=nothing, tmp_dir=nothing, render_svg=true)

Render the standard solution form from the final augmented matrix.
"""
function show_solution( matrices; var_name::String="x", fig_scale=1, output_dir=nothing, tmp_dir=nothing, render_svg=true, render_opts=nothing )
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
    if render_svg
        return _render_solution_svg(tex; fig_scale=fig_scale, output_dir=_resolve_output_dir(output_dir, tmp_dir), render_opts=render_opts)
    end
    return _display_tex(tex)
end
