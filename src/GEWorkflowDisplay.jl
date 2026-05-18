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
    la = load_LAFigureSpecs()
    linear_system_tex = _pygetattr(la, :linear_system_tex)
    tex = _pycall(linear_system_tex, A, b; var_name=var_name)
    py = _ensure_pythoncall()
    tex = Base.invokelatest(py.pyconvert, String, tex)
    bs = _pyimport("matrixlayout.backsubst")
    backsubst_svg = _pygetattr(bs, :backsubst_svg)
    resolved_output_dir = _resolve_output_dir(output_dir, tmp_dir, pb.tmp_dir)
    svg = _pycall(backsubst_svg; system_txt=tex, show_system=true,
                  show_cascade=false, show_solution=false,
                  fig_scale=fig_scale, output_dir=resolved_output_dir,
                  render_opts=render_opts)
    svg_str = Base.invokelatest(py.pyconvert, String, svg)
    return SVGOut(svg_str)
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
        py = _ensure_pythoncall()
        lines = Base.invokelatest(py.pyconvert, Vector{String}, lines)
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
        py = _ensure_pythoncall()
        lines = Base.invokelatest(py.pyconvert, Vector{String}, lines)
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
        py = _ensure_pythoncall()
        tex = Base.invokelatest(py.pyconvert, String, tex)
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
    lines = load_LAFigureSpecs().backsubstitution_tex(A, b, var_name=var_name)
    return _render_backsubst_svg(lines; fig_scale=fig_scale, output_dir=pb.tmp_dir, render_opts=render_opts)
end
# --------------------------------------------------------------------------------------------------------------
function show_forwardsubstitution!(  pb::ShowGE{T}; b_col=1, var_name::String="x", fig_scale=1, render_svg=true, render_opts=nothing ) where T <: Number
    A, b = _forwardsub_ref(pb; b_col=b_col)
    lines = load_LAFigureSpecs().backsubstitution_tex(A[end:-1:1, end:-1:1], b[end:-1:1], var_name=var_name)
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
    tex = load_LAFigureSpecs().standard_solution_tex(A, b, var_name=var_name)
    return _render_solution_svg(tex; fig_scale=fig_scale, render_opts=render_opts)
end
# ==============================================================================================================
raw"""
    show_backsubstitution(A, b; var_name="x", fig_scale=1, output_dir=nothing, tmp_dir="/tmp/la/run", render_svg=true)

    Render the back-substitution cascade for the upper-triangular system `A * x = b`
    using `LAFigureSpecs.backsubstitution_tex`. Works with Integer/Float as well as
    exact `Rational` and `Complex{Rational}` inputs (those are converted to tuples so
    SymPy reconstructs exact rationals on the Python side).
"""
function show_backsubstitution(A, b; var_name::String="x", fig_scale=1, output_dir=nothing, tmp_dir="/tmp/la/run", render_svg=true, render_opts=nothing)
    A2 = (A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(A) : A
    b2 = (b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(b) : b
    lines = load_LAFigureSpecs().backsubstitution_tex(A2, b2, var_name=var_name)
    if render_svg
        return _render_backsubst_svg(lines; fig_scale=fig_scale, output_dir=_resolve_output_dir(output_dir, tmp_dir), render_opts=render_opts)
    end
    return _display_cascade(lines)
end
# --------------------------------------------------------------------------------------------------------------
raw"""
    show_forwardsubstitution(A, b; var_name="x", fig_scale=1, output_dir=nothing, tmp_dir="/tmp/la/run", render_svg=true)

Render the forward-substitution cascade for the lower-triangular system `A * x = b`
using the LAFigureSpecs backsubstitution cascade on a reversed system, then relabeling indices.
Supports Integer/Float as well as exact `Rational` and `Complex{Rational}` inputs
converted to tuples for exact SymPy reconstruction.
"""
function show_forwardsubstitution(A, b; var_name::String="x", fig_scale=1, output_dir=nothing, tmp_dir="/tmp/la/run", render_svg=true, render_opts=nothing)
    A2 = (A isa AbstractArray{<:Rational} || A isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(A) : A
    b2 = (b isa AbstractArray{<:Rational} || b isa AbstractArray{Complex{<:Rational}}) ? _encode_exact.(b) : b
    lines = load_LAFigureSpecs().backsubstitution_tex(A2[end:-1:1, end:-1:1], b2[end:-1:1], var_name=var_name)
    lines = _relabel_cascade(lines, size(A, 1); var_name=var_name)
    if render_svg
        return _render_backsubst_svg(lines; fig_scale=fig_scale, output_dir=_resolve_output_dir(output_dir, tmp_dir), render_opts=render_opts)
    end
    return _display_cascade(lines)
end
# ==============================================================================================================
