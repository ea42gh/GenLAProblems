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
