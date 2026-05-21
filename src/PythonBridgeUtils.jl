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

    for T in (String, Int, Float64, Bool)
        try
            return Base.invokelatest(py.pyconvert, T, x)
        catch
        end
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
        tup = Base.invokelatest(py.pyconvert, Tuple, x)
        return tuple((materialize_python_value(v) for v in tup)...)
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
function _is_missing_python_module(err, module_name::AbstractString)
    msg = sprint(showerror, err)
    return occursin("ModuleNotFoundError", msg) && occursin("No module named '$module_name'", msg)
end

function load_LAFigureSpecs()
    if _LAFigureSpecs[] === nothing
        pc = _ensure_pythoncall()
        if pc === nothing
            return nothing
        end
        try
            _LAFigureSpecs[] = Base.invokelatest(pc.pyimport, "LAFigureSpecs")
        catch err
            if _is_missing_python_module(err, "LAFigureSpecs")
                error(
                    "Python module `LAFigureSpecs` is required by GenLAProblems.\n" *
                    "Install it in the active Python environment.\n\n" *
                    "Original error:\n$err"
                )
            end
            rethrow(err)
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
        pc = _ensure_pythoncall()
        if pc === nothing
            return nothing
        end
        try
            _matrixlayout[] = Base.invokelatest(pc.pyimport, "matrixlayout")
        catch err
            if _is_missing_python_module(err, "matrixlayout")
                error(
                    "Python module `matrixlayout` is required by GenLAProblems.\n" *
                    "Install it in the active Python environment.\n\n" *
                    "Original error:\n$err"
                )
            end
            rethrow(err)
        end
    end
    return _matrixlayout[]
end
