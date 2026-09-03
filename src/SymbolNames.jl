# ------------------------------------------------------------------------------
# ---------------------------------------------- matrices and vectors of symbols
# ------------------------------------------------------------------------------
"""
    symbol_vector(s, indices) -> Vector{Symbol}

Construct a vector of symbolic names using the base string `s` and the supplied
indices.

For example, `symbol_vector("x", 1:3)` returns `[:x_1, :x_2, :x_3]`.
"""
function _symbol_base_string(s)
    sstr = string(s)
    if startswith(sstr, "\$") && endswith(sstr, "\$") && length(sstr) >= 2
        return sstr[2:(end-1)]
    end
    return sstr
end

function symbol_vector(s, indices)
    sstr = _symbol_base_string(s)
    [Symbol(sstr * "_$i") for i in collect(indices)]
end
# ------------------------------------------------------------------------------
"""
    symbols_matrix(s, row_indices, col_indices) -> Matrix{Symbol}

Create a matrix of symbolic names using the base string `s`.
Entry `(i, j)` is named like `:a_{1,3}` after collecting the supplied row and
column index iterables.

For example, `symbols_matrix("a", 1:2, 3:4)` returns
`[:a_{1,3} :a_{1,4}; :a_{2,3} :a_{2,4}]`.
"""
function symbols_matrix(s, row_indices, col_indices)
    sstr = _symbol_base_string(s)
    rows = collect(row_indices)
    cols = collect(col_indices)
    matrix = [Symbol("$(sstr)_{$(i),$(j)}") for i in rows, j in cols]
    return matrix
end
