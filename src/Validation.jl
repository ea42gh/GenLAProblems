# ------------------------------------------------------------------------------
# Shared validation and integer sampling helpers
# ------------------------------------------------------------------------------
"""
    _int_range(maxint, has_zeros) -> AbstractVector{Int}

Return the integer range used to populate random matrices, optionally including zero.
"""
function _validate_dimension(fname, name, value; minimum = 0)
    value isa Integer && value >= minimum ||
        throw(ArgumentError("$fname requires $name to be an integer >= $minimum"))
end

function _validate_block_sizes(fname, blocks)
    values = blocks isa Integer ? [blocks] : collect(blocks)
    (!isempty(values) && all(x -> x isa Integer && x > 0, values)) || throw(
        ArgumentError(
            "$fname requires a nonempty collection of positive integer block sizes",
        ),
    )
    values
end

function _validate_maxint(fname, maxint)
    maxint isa Integer && maxint >= 0 ||
        throw(ArgumentError("$fname requires maxint to be a nonnegative integer"))
end
function _validate_positive_maxint(fname, maxint)
    _validate_maxint(fname, maxint)
    maxint > 0 || throw(ArgumentError("$fname requires maxint to be positive"))
end

function _int_range(maxint, has_zeros)
    _validate_maxint("integer matrix generator", maxint)
    if has_zeros
        rng = (-maxint):maxint
    else
        rng = [(-maxint):-1; 1:maxint]
    end
    rng
end

function _validate_rank_request(fname, m, n, r)
    all(x -> x isa Integer, (m, n, r)) ||
        throw(ArgumentError("$fname requires integer dimensions and rank"))
    m >= 0 && n >= 0 ||
        throw(ArgumentError("$fname requires nonnegative matrix dimensions"))
    0 <= r <= min(m, n) ||
        throw(ArgumentError("$fname requires r to satisfy 0 <= r <= min(m, n)"))
end
