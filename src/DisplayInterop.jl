"""
    SVGOut(svg::String)

Wrapper type for SVG output used by notebook-friendly display helpers.
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

function _named_qr_matrices(qr)
    getmat(name::String) = begin
        if qr isa NamedTuple
            return getproperty(qr, Symbol(name))
        elseif qr isa AbstractDict
            return get(qr, name, get(qr, Symbol(name), nothing))
        end
        return nothing
    end

    return (
        A = getmat("A"),
        W = getmat("W"),
        WtA = getmat("WtA"),
        WtW = getmat("WtW"),
        S = getmat("S"),
        Qt = getmat("Qt"),
        Q = getmat("Q"),
        R = getmat("R"),
    )
end

function _show_svg(svg)
    if svg === nothing
        return SVGOut("")
    end
    if _py_is_py(svg)
        _ensure_pythoncall()
        py = getfield(@__MODULE__, :PythonCall)
        if _py_is_none(svg)
            return SVGOut("")
        end
        return SVGOut(Base.invokelatest(py.pyconvert, String, svg))
    elseif svg isa AbstractString
        return SVGOut(svg)
    elseif isdefined(@__MODULE__, :PythonCall)
        py = getfield(@__MODULE__, :PythonCall)
        try
            return SVGOut(Base.invokelatest(py.pyconvert, String, svg))
        catch
        end
    end
    return svg
end

"""
    py_show_svg(svg)

Display an SVG in a Python notebook (e.g., %%julia cell) via IPython.display.SVG.
Accepts `SVGOut`, raw SVG strings, or PythonCall.Py SVG objects.
"""
function py_show_svg(svg)
    try
        _ensure_pythoncall()
        ip = _pyimport("IPython.display")
        py_display = _pygetattr(ip, :display)
        py_svg = _pygetattr(ip, :SVG)
        if svg isa SVGOut
            return _pycall(py_display, _pycall(py_svg, svg.svg))
        elseif svg isa AbstractString
            return _pycall(py_display, _pycall(py_svg, svg))
        elseif _py_is_py(svg)
            s = Base.invokelatest(getfield(@__MODULE__, :PythonCall).pyconvert, String, svg)
            return _pycall(py_display, _pycall(py_svg, s))
        else
            error("py_show_svg expects SVGOut, SVG string, or PythonCall.Py")
        end
    catch
        if svg isa SVGOut
            return Base.display(svg)
        elseif svg isa AbstractString
            return Base.display(SVGOut(svg))
        end
        error("py_show_svg expects SVGOut, SVG string, or PythonCall.Py")
    end
end

"""
    show_svg(svg)

Alias for `py_show_svg`, for notebook-friendly SVG display.
"""
show_svg(svg) = py_show_svg(svg)

"""
    l_show_svd(A, U, Σ, Vt, rankA)

Display an SVD factorization with block structure separating the rank and null
space components.
"""
function l_show_svd(A, U, Σ, Vt, rankA)
    BA = _ensure_blockarrays()
    Ub = BA.BlockArray(U, [size(U, 1)], [rankA, size(U, 2) - rankA])
    Σb = BA.BlockArray(Σ, [rankA, size(Σ, 1) - rankA], [rankA, size(Σ, 2) - rankA])
    Vtb = BA.BlockArray(Vt, [rankA, size(Vt, 1) - rankA], [size(Vt, 2)])
    display(l_show(LaTeXString("A = U \\Sigma V^T : "), A, " = ", Ub, Σb, Vtb))
    return nothing
end

function _clean_tmp_kwargs(kwargs)
    clean = Dict(kwargs)
    pop!(clean, :tmp_dir, nothing)
    pop!(clean, :keep_file, nothing)
    pop!(clean, :output_dir, nothing)
    return clean
end

function _map_tmp_to_output(kwargs)
    clean = Dict(kwargs)
    if haskey(clean, :tmp_dir) && !haskey(clean, :output_dir)
        clean[:output_dir] = clean[:tmp_dir]
    end
    pop!(clean, :tmp_dir, nothing)
    pop!(clean, :keep_file, nothing)
    return clean
end

function _normalize_render_opts(render_opts; tmp_dir=nothing)
    opts = if render_opts === nothing
        Dict{String,Any}()
    elseif !(render_opts isa AbstractDict)
        Dict{String,Any}(render_opts)
    else
        Dict{String,Any}(render_opts)
    end
    if tmp_dir !== nothing && !haskey(opts, "output_dir")
        opts["output_dir"] = tmp_dir
    end
    return opts
end
