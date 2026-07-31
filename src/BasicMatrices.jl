# ------------------------------------------------------------------------------
# ------------------------------------- matrices for use in GE and GJ algorithms
# ------------------------------------------------------------------------------
"""
    _int_range(maxint, has_zeros) -> AbstractVector{Int}

Return the integer range used to populate random matrices, optionally including zero.
"""
function _int_range( maxint, has_zeros)
    if has_zeros
        rng = -maxint:maxint
    else
        rng = [-maxint:-1; 1:maxint]
    end
    rng
end

function _validate_rank_request(fname, m, n, r)
    0 <= r <= min(m, n) || throw(ArgumentError("$fname requires r to satisfy 0 <= r <= min(m, n)"))
end
# ------------------------------------------------------------------------------
"""
    unit_lower(m, n; maxint=3) -> Matrix{Int}
    unit_lower(m; maxint=3) -> Matrix{Int}

Construct an `m × n` unit lower-triangular integer matrix.
Entries strictly below the main diagonal are sampled from `-maxint:maxint`,
diagonal entries are `1`, and entries above the diagonal are `0`.
The one-argument method returns the square `m × m` case.
"""
function unit_lower(m,n; maxint=3)
    # create a unit lower triangular matrix
    [ x>y ? rand(-maxint:maxint) : (x == y ? 1 : 0) for x in 1:m, y in 1:n]
end
# ------------------------------------------------------------------------------
function unit_lower(m; maxint=3)
   unit_lower(m,m,maxint=maxint)
end
# ------------------------------------------------------------------------------
"""
    lower(m, n; maxint=3) -> Matrix{Int}
    lower(m; maxint=3) -> Matrix{Int}

Construct an integer lower-triangular matrix with nonzero diagonal.
This starts from `unit_lower` and replaces each diagonal entry by a random
integer chosen from `[-maxint:-1; 1:maxint]`.
The one-argument method returns the square `m × m` case.
"""
function lower(m,n; maxint=3)
    L = unit_lower(m,n; maxint=maxint)
    for i in 1:min(m,n)
        L[i,i] = rand( [-maxint:-1; 1:maxint])
    end
    L
end
# ------------------------------------------------------------------------------
function lower(m; maxint=3)
    lower(m,m,maxint=maxint)
end
# ------------------------------------------------------------------------------
"""
    rref_matrix(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false) -> R, pivot_cols

Generate an `m × n` reduced row-echelon matrix of rank `r`.
`pivot_cols` contains the pivot column indices. When `pivot_in_first_col=true`,
column `1` is forced to be a pivot unless `r == 0`. Nonpivot entries are
sampled from a small integer range controlled by `maxint`; zeros may also be
sampled when `has_zeros=true`.
"""
function rref_matrix(m,n,r; maxint=3, pivot_in_first_col=true, has_zeros=false)
    # create a reduced row echelon form matrix of size m x n and rank r
    _validate_rank_request("rref_matrix", m, n, r)
    if pivot_in_first_col || r==n
        pivot_cols = sort!([1; (2:n)[randperm(n-1)]][1:r])
    else
        pivot_cols = sort!((2:n)[randperm(n-1)][1:r])
    end

    rng = _int_range(maxint,has_zeros)

    if m > r
        M = [ rand(rng, (r,n))
              zeros(Int64, (m-r,n))
        ]
    else
        M = rand( rng, (m,n) )
    end
    for i in 1:r
        for j in 1:(pivot_cols[i]-1)
            M[i,j] = 0
        end
        M[i,pivot_cols[i]]         = 1
        M[1:(i-1), pivot_cols[i]] .= 0
    end
    M, pivot_cols
end
# ------------------------------------------------------------------------------
"""
    ref_matrix(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false) -> U, pivot_cols

Generate an `m × n` row-echelon-form matrix of rank `r`.
This starts from `rref_matrix(...)` and then scales pivot rows and mixes
nonpivot columns so the result remains in row echelon form rather than reduced
row echelon form. `pivot_cols` lists the pivot columns.
"""
function ref_matrix(m,n,r; maxint=3, pivot_in_first_col=true, has_zeros=false)
    _validate_rank_request("ref_matrix", m, n, r)
    M,pivot_cols = rref_matrix(m,n,r; maxint=maxint, pivot_in_first_col=pivot_in_first_col, has_zeros=has_zeros)
    rng = _int_range( maxint, false)
    M   = Diagonal(rand(rng,m)) * M * unit_lower(n,n,maxint=1)'
    M, pivot_cols
end
# ------------------------------------------------------------------------------
"""
    gen_full_col_rank_matrix(mc, nc; maxint=3) -> Matrix{Int}

Generate a small exact matrix with full column rank.
`mc` and `nc` may be integers or block-size collections whose sums define the
row and column counts. The construction uses a structured orthogonal factor and
an invertible lower-triangular mixing matrix so that the resulting matrix has
rank equal to its number of columns.
"""
function gen_full_col_rank_matrix(mc,nc; maxint=3)
    # produce a reasonable A'A matrix; need m ≥ n
    m = sum(mc)
    n = sum(nc)

    Q = sparse_Q_matrix(mc)
    M = zeros(Int64, (m,n))
    rng = _int_range(maxint, false)
    for i = 1:min(m,n)
        M[i,i] = rand( rng )
    end
    Q[:,1:m]*unit_lower(m,maxint=maxint)*M
end
# ------------------------------------------------------------------------------
"""
    symmetric_matrix(m; maxint=3, with_zeros=false) -> Matrix{Int}

Generate a random symmetric `m × m` integer matrix with nonzero diagonal.
Off-diagonal entries are sampled from the integer range controlled by `maxint`;
when `with_zeros=true`, zeros may also appear off the diagonal.
"""
function symmetric_matrix(m;maxint=3, with_zeros=false )
    rng = _int_range(maxint,with_zeros)
    A = [ x>y ? rand(rng) : 0 for x in 1:m, y in 1:m]
    A = A+A'
    for i in 1:m
        A[i,i] = rand([-maxint:-1; 1:maxint])
    end
    A
end
# ------------------------------------------------------------------------------
"""
    skew_symmetric_matrix(m; maxint=5, with_zeros=false) -> Matrix{Int}

Generate a random skew-symmetric `m × m` integer matrix.
The diagonal is zero by construction, and entries below the diagonal are drawn
from the integer range controlled by `maxint`.
"""
function skew_symmetric_matrix(m;maxint=5, with_zeros=false )
    rng = _int_range(maxint,with_zeros)
    A = [ i>j ? rand(rng) : 0 for i in 1:m, j in 1:m]
    A - A'
end
# ------------------------------------------------------------------------------
"""
    e_i(i, n) -> Vector{Int}

Return the `i`th standard basis vector in `R^n`.
"""
function e_i(i,n)
    v = zeros( Int, n )
    v[i] = 1
    v
end
# ------------------------------------------------------------------------------
"""
    i_with_onecol(m, c; maxint=3, with_zeros=false, lower=true, upper=true) -> Matrix{Int}

Construct an identity-like elimination matrix whose `c`th column is randomized.
The diagonal entry `(c,c)` remains `1`. Entries below that position are filled
when `lower=true`, entries above it are filled when `upper=true`, and sampled
values are controlled by `maxint` and `with_zeros`.
"""
function i_with_onecol(m,c; maxint=3, with_zeros=false, lower=true, upper=true)
    rng = _int_range(maxint,with_zeros)
    # take I and set column c to random entries
    E        = collect(1I(m))           # Int64  eye(m)
    if lower && c < m
        E[c+1:m,c] = rand(rng, m-c)  # set column c to non-zero entries
    end
    if upper && c > 1
        E[1:c-1,c] = rand(rng, c-1)  # set column c to non-zero entries
    end
    E[c,c]   = 1
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
    sort(collect(row_order)) == collect(1:n) || throw(ArgumentError("row_order must be a permutation of 1:$n"))
    row_order == collect(1:n) && throw(ArgumentError("row_order must not be the identity permutation"))
    P = zeros(Int, (n,n))
    for i in 1:n
        P[row_order[i],i] = 1
    end
    P
end
# ------------------------------------------------------------------------------
"""
    gen_permutation_matrix(n) -> Matrix{Int}

Construct a uniformly shuffled `n × n` permutation matrix.
"""
function gen_permutation_matrix(n)
    n > 1 || throw(ArgumentError("n must be at least 2 to generate a non-identity permutation"))
    locs = randperm(n)
    while locs == collect(1:n)
        locs = randperm(n)
    end
    P    = zeros(Int, (n,n))
    for i in 1:n
        P[i,locs[i]] = 1
    end
    P
end
