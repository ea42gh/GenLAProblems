# ------------------------------------------------------------------------------
# Basis and permutation matrix utilities
# ------------------------------------------------------------------------------
"""
    e_i(i, n) -> Vector{Int}

Return the `i`th standard basis vector in `R^n`.
"""
function e_i(i, n)
    _validate_dimension("e_i", "n", n; minimum = 1)
    1 <= i <= n || throw(ArgumentError("e_i requires i to satisfy 1 <= i <= n"))
    v = zeros(Int, n)
    v[i] = 1
    v
end
# ------------------------------------------------------------------------------
"""
    i_with_onecol(m, c; maxint=3, with_zeros=false, lower=true, upper=true, rng=Random.default_rng()) -> Matrix{Int}

Construct an identity-like elimination matrix whose `c`th column is randomized.
The diagonal entry `(c,c)` remains `1`. Entries below that position are filled
when `lower=true`, entries above it are filled when `upper=true`, and sampled
values are controlled by `maxint` and `with_zeros`.
"""
function i_with_onecol(
    m,
    c;
    maxint = 3,
    with_zeros = false,
    lower = true,
    upper = true,
    rng = Random.default_rng(),
)
    _validate_dimension("i_with_onecol", "m", m; minimum = 1)
    1 <= c <= m || throw(ArgumentError("i_with_onecol requires c to satisfy 1 <= c <= m"))
    _validate_maxint("i_with_onecol", maxint)
    needs_samples = (lower && c < m) || (upper && c > 1)
    values = needs_samples ? _int_range(maxint, with_zeros) : Int[]
    needs_samples && !with_zeros && _validate_positive_maxint("i_with_onecol", maxint)
    # take I and set column c to random entries
    E = collect(1I(m))           # Int64  eye(m)
    if lower && c < m
        E[(c+1):m, c] = rand(rng, values, m - c)  # set column c to non-zero entries
    end
    if upper && c > 1
        E[1:(c-1), c] = rand(rng, values, c - 1)  # set column c to non-zero entries
    end
    E[c, c] = 1
    E
end
# ------------------------------------------------------------------------------
"""
    gen_permutation_matrix(row_order::AbstractVector{<:Integer}) -> Matrix{Int}

Construct the permutation matrix determined by `row_order`.
If `row_order[i] == j`, then column `i` of the result is the basis vector `e_j`.
"""
function gen_permutation_matrix(row_order::AbstractVector{<:Integer})
    n = length(row_order)
    sort(collect(row_order)) == collect(1:n) ||
        throw(ArgumentError("row_order must be a permutation of 1:$n"))
    row_order == collect(1:n) &&
        throw(ArgumentError("row_order must not be the identity permutation"))
    P = zeros(Int, (n, n))
    for i = 1:n
        P[row_order[i], i] = 1
    end
    P
end
# ------------------------------------------------------------------------------
"""
    gen_permutation_matrix(n; rng=Random.default_rng()) -> Matrix{Int}

Construct a uniformly shuffled `n × n` permutation matrix.
"""
function gen_permutation_matrix(n; rng = Random.default_rng())
    _validate_dimension("gen_permutation_matrix", "n", n; minimum = 2)
    n > 1 ||
        throw(ArgumentError("n must be at least 2 to generate a non-identity permutation"))
    locs = randperm(rng, n)
    while locs == collect(1:n)
        locs = randperm(rng, n)
    end
    P = zeros(Int, (n, n))
    for i = 1:n
        P[i, locs[i]] = 1
    end
    P
end
