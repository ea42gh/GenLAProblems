"""
    la_version(), la_build(), ml_version(), ml_build()

Return version/build strings exposed by the Python packages `LAFigureSpecs` and
`matrixlayout`. These require PythonCall at runtime.
"""
function la_version()
    py = _ensure_pythoncall()
    la = load_LAFigureSpecs()
    v = _pygetattr(la, :__version__)
    return Base.invokelatest(py.pyconvert, String, v)
end

function la_build()
    py = _ensure_pythoncall()
    la = load_LAFigureSpecs()
    v = _pygetattr(la, :__build__)
    return Base.invokelatest(py.pyconvert, String, v)
end

function ml_version()
    py = _ensure_pythoncall()
    ml = load_matrixlayout()
    v = _pygetattr(ml, :__version__)
    return Base.invokelatest(py.pyconvert, String, v)
end

function ml_build()
    py = _ensure_pythoncall()
    ml = load_matrixlayout()
    v = _pygetattr(ml, :__build__)
    return Base.invokelatest(py.pyconvert, String, v)
end

struct SympyProxy end
const sympy = SympyProxy()

struct NMProxy end
const nM = NMProxy()

(::NMProxy)() = nM

function materialize_python_value(x)
    if x isa NamedTuple
        return (; (name => materialize_python_value(value) for (name, value) in pairs(x))...)
    elseif !_py_is_py(x)
        if x isa AbstractDict
            return Dict(materialize_python_value(k) => materialize_python_value(v) for (k, v) in pairs(x))
        elseif x isa Tuple
            return tuple((materialize_python_value(v) for v in x)...)
        elseif x isa AbstractArray
            return map(materialize_python_value, x)
        end
        return x
    end

    if _py_is_none(x)
        return nothing
    end

    py = _ensure_pythoncall()
    converted = try
        Base.invokelatest(py.pyconvert, Any, x)
    catch
        x
    end

    if converted !== x
        return materialize_python_value(converted)
    end

    try
        shape = Base.invokelatest(py.pygetattr, x, "shape")
        shp = Base.invokelatest(py.pyconvert, Tuple, shape)
        if length(shp) == 2
            return map(materialize_python_value, sym_to_julia_mat(x))
        elseif length(shp) == 1
            return map(materialize_python_value, sym_to_julia_vec(x))
        end
    catch
    end

    try
        items_fn = Base.invokelatest(py.pygetattr, x, "items")
        items = _pycall(items_fn)
        pairs_vec = Base.invokelatest(py.pyconvert, Vector{Any}, items)
        return Dict(
            begin
                pair_t = Base.invokelatest(py.pyconvert, Tuple, pair)
                materialize_python_value(pair_t[1]) => materialize_python_value(pair_t[2])
            end for pair in pairs_vec
        )
    catch
    end

    try
        vec = Base.invokelatest(py.pyconvert, Vector{Any}, x)
        return map(materialize_python_value, vec)
    catch
    end

    return x
end

"""
    load_LAFigureSpecs() -> LAFigureSpecs

Load the Python `LAFigureSpecs` module via PythonCall.
"""
function load_LAFigureSpecs()
    if _LAFigureSpecs[] === nothing
        try
            pc = _ensure_pythoncall()
            if pc === nothing
                return nothing
            end
            _LAFigureSpecs[] = Base.invokelatest(pc.pyimport, "LAFigureSpecs")
        catch err
            error(
                "Python module `LAFigureSpecs` is required by GenLAProblems.\n" *
                "Install it in the active Python environment.\n\n" *
                "Original error:\n$err"
            )
        end
    end
    return _LAFigureSpecs[]
end

"""
    load_matrixlayout() -> matrixlayout

Load the Python `matrixlayout` module via PythonCall.
"""
function load_matrixlayout()
    if _matrixlayout[] === nothing
        try
            pc = _ensure_pythoncall()
            if pc === nothing
                return nothing
            end
            _matrixlayout[] = Base.invokelatest(pc.pyimport, "matrixlayout")
        catch err
            error(
                "Python module `matrixlayout` is required by GenLAProblems.\n" *
                "Install it in the active Python environment.\n\n" *
                "Original error:\n$err"
            )
        end
    end
    return _matrixlayout[]
end

"""
    nM -> NMProxy

Proxy that exposes matrixlayout helpers by default and LAFigureSpecs display helpers.
"""
