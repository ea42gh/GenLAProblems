# ------------------------------------------------------------------------------
# ------------------------------------- matrices for use in GE and GJ algorithms
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

"""
    gen_matrix(m, n; maxint=3, minint=-maxint, rng=Random.default_rng()) -> Matrix{Int}

Construct an unconstrained `m × n` integer matrix by sampling each entry
uniformly from `minint:maxint`. This generator does not impose a rank or
other structural condition.
"""
function gen_matrix(m, n; maxint = 3, minint = -maxint, rng = Random.default_rng())
    _validate_dimension("gen_matrix", "m", m)
    _validate_dimension("gen_matrix", "n", n)
    _validate_maxint("gen_matrix", maxint)
    minint isa Integer && minint <= maxint ||
        throw(ArgumentError("gen_matrix requires minint to be an integer <= maxint"))
    rand(rng, minint:maxint, m, n)
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
# ------------------------------------------------------------------------------
"""
    unit_lower(m, n; maxint=3, rng=Random.default_rng()) -> Matrix{Int}
    unit_lower(m; maxint=3, rng=Random.default_rng()) -> Matrix{Int}

Construct an `m × n` unit lower-triangular integer matrix.
Entries strictly below the main diagonal are sampled from `-maxint:maxint`,
diagonal entries are `1`, and entries above the diagonal are `0`.
The one-argument method returns the square `m × m` case.
"""
function unit_lower(m, n; maxint = 3, rng = Random.default_rng())
    _validate_dimension("unit_lower", "m", m)
    _validate_dimension("unit_lower", "n", n)
    _validate_maxint("unit_lower", maxint)
    # create a unit lower triangular matrix
    [x > y ? rand(rng, (-maxint):maxint) : (x == y ? 1 : 0) for x = 1:m, y = 1:n]
end
# ------------------------------------------------------------------------------
function unit_lower(m; maxint = 3, rng = Random.default_rng())
    unit_lower(m, m, maxint = maxint, rng = rng)
end
# ------------------------------------------------------------------------------
"""
    lower(m, n; maxint=3, rng=Random.default_rng()) -> Matrix{Int}
    lower(m; maxint=3, rng=Random.default_rng()) -> Matrix{Int}

Construct an integer lower-triangular matrix with nonzero diagonal.
This starts from `unit_lower` and replaces each diagonal entry by a random
integer chosen from `[-maxint:-1; 1:maxint]`.
The one-argument method returns the square `m × m` case.
"""
function lower(m, n; maxint = 3, rng = Random.default_rng())
    _validate_maxint("lower", maxint)
    min(m, n) == 0 || _validate_positive_maxint("lower", maxint)
    L = unit_lower(m, n; maxint = maxint, rng = rng)
    for i = 1:min(m, n)
        L[i, i] = rand(rng, [(-maxint):-1; 1:maxint])
    end
    L
end
# ------------------------------------------------------------------------------
function lower(m; maxint = 3, rng = Random.default_rng())
    lower(m, m, maxint = maxint, rng = rng)
end
# ------------------------------------------------------------------------------
"""
    rref_matrix(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false, rng=Random.default_rng()) -> R, pivot_cols

Generate an `m × n` reduced row-echelon matrix of rank `r`.
`pivot_cols` contains the pivot column indices. When `pivot_in_first_col=true`,
column `1` is forced to be a pivot unless `r == 0`. Nonpivot entries are
sampled from a small integer range controlled by `maxint`; zeros may also be
sampled when `has_zeros=true`.
"""
function rref_matrix(
    m,
    n,
    r;
    maxint = 3,
    pivot_in_first_col = true,
    has_zeros = false,
    rng = Random.default_rng(),
)
    # create a reduced row echelon form matrix of size m x n and rank r
    _validate_rank_request("rref_matrix", m, n, r)
    _validate_maxint("rref_matrix", maxint)
    r == 0 && return zeros(Int64, m, n), Int[]
    (has_zeros || r == 0) || _validate_positive_maxint("rref_matrix", maxint)
    if pivot_in_first_col || r == n
        pivot_cols = sort!([1; (2:n)[randperm(rng, n - 1)]][1:r])
    else
        pivot_cols = sort!((2:n)[randperm(rng, n - 1)][1:r])
    end

    values = _int_range(maxint, has_zeros)

    if m > r
        M = [
            rand(rng, values, (r, n))
            zeros(Int64, (m - r, n))
        ]
    else
        M = rand(rng, values, (m, n))
    end
    for i = 1:r
        for j = 1:(pivot_cols[i]-1)
            M[i, j] = 0
        end
        M[i, pivot_cols[i]] = 1
        M[1:(i-1), pivot_cols[i]] .= 0
    end
    M, pivot_cols
end
# ------------------------------------------------------------------------------
"""
    ref_matrix(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false, rng=Random.default_rng()) -> U, pivot_cols

Generate an `m × n` row-echelon-form matrix of rank `r`.
This starts from `rref_matrix(...)` and then scales pivot rows and mixes
nonpivot columns so the result remains in row echelon form rather than reduced
row echelon form. `pivot_cols` lists the pivot columns.
"""
function ref_matrix(
    m,
    n,
    r;
    maxint = 3,
    pivot_in_first_col = true,
    has_zeros = false,
    rng = Random.default_rng(),
)
    _validate_rank_request("ref_matrix", m, n, r)
    _validate_maxint("ref_matrix", maxint)
    r == 0 && return zeros(Int64, m, n), Int[]
    m == 0 || _validate_positive_maxint("ref_matrix", maxint)
    M, pivot_cols = rref_matrix(
        m,
        n,
        r;
        maxint = maxint,
        pivot_in_first_col = pivot_in_first_col,
        has_zeros = has_zeros,
        rng = rng,
    )
    values = _int_range(maxint, false)
    M = Diagonal(rand(rng, values, m)) * M * unit_lower(n, n, maxint = 1, rng = rng)'
    M, pivot_cols
end
# ------------------------------------------------------------------------------
"""
    gen_full_col_rank_matrix(mc, nc; maxint=3, rng=Random.default_rng()) -> Matrix{Int}

Generate a small exact matrix with full column rank.
`mc` and `nc` may be integers or block-size collections whose sums define the
row and column counts. The construction uses a structured orthogonal factor and
an invertible lower-triangular mixing matrix so that the resulting matrix has
rank equal to its number of columns.
"""
function gen_full_col_rank_matrix(mc, nc; maxint = 3, rng = Random.default_rng())
    # produce a reasonable A'A matrix; need m ≥ n
    row_blocks = _validate_block_sizes("gen_full_col_rank_matrix", mc)
    col_blocks = _validate_block_sizes("gen_full_col_rank_matrix", nc)
    m = sum(row_blocks)
    n = sum(col_blocks)
    m >= n || throw(ArgumentError("gen_full_col_rank_matrix requires m >= n"))
    _validate_positive_maxint("gen_full_col_rank_matrix", maxint)

    Q = sparse_Q_matrix(row_blocks; rng = rng)
    M = zeros(Int64, (m, n))
    values = _int_range(maxint, false)
    for i = 1:min(m, n)
        M[i, i] = rand(rng, values)
    end
    Q[:, 1:m] * unit_lower(m, maxint = maxint, rng = rng) * M
end
# ------------------------------------------------------------------------------
"""
    symmetric_matrix(m; maxint=3, with_zeros=false, rng=Random.default_rng()) -> Matrix{Int}

Generate a random symmetric `m × m` integer matrix with nonzero diagonal.
Off-diagonal entries are sampled from the integer range controlled by `maxint`;
when `with_zeros=true`, zeros may also appear off the diagonal.
"""
function symmetric_matrix(m; maxint = 3, with_zeros = false, rng = Random.default_rng())
    _validate_dimension("symmetric_matrix", "m", m)
    _validate_maxint("symmetric_matrix", maxint)
    m == 0 || _validate_positive_maxint("symmetric_matrix", maxint)
    values = _int_range(maxint, with_zeros)
    A = [x > y ? rand(rng, values) : 0 for x = 1:m, y = 1:m]
    A = A + A'
    for i = 1:m
        A[i, i] = rand(rng, [(-maxint):-1; 1:maxint])
    end
    A
end
# ------------------------------------------------------------------------------
"""
    skew_symmetric_matrix(m; maxint=5, with_zeros=false, rng=Random.default_rng()) -> Matrix{Int}

Generate a random skew-symmetric `m × m` integer matrix.
The diagonal is zero by construction, and entries below the diagonal are drawn
from the integer range controlled by `maxint`.
"""
function skew_symmetric_matrix(
    m;
    maxint = 5,
    with_zeros = false,
    rng = Random.default_rng(),
)
    _validate_dimension("skew_symmetric_matrix", "m", m)
    _validate_maxint("skew_symmetric_matrix", maxint)
    (with_zeros || m <= 1) || _validate_positive_maxint("skew_symmetric_matrix", maxint)
    values = _int_range(maxint, with_zeros)
    A = [i > j ? rand(rng, values) : 0 for i = 1:m, j = 1:m]
    A - A'
end
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
function gen_permutation_matrix(
    row_order::AbstractVector{<:Integer};
    rng = Random.default_rng(),
)
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
