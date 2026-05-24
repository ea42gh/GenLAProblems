#using LinearAlgebra, Latexify, Random

# ------------------------------------------------------------------------------
# ----------- Integer Square Roots:  Pythagorean Number Triplets and Quadruplets
# ------------------------------------------------------------------------------
PythagoreanNumberTriplets =
[   3    4    5
    5   12   13
    7   24   25
    8   15   17
    9   40   41
   11   60   61
   12   35   37
   13   84   85
   15  112  113 ]

PythagoreanNumberQuadruplets =
[   1   2   2   3
    2  10  11  15
    4  13  16  21
    2  10  25  27
    2   3   6   7
    1  12  12  17
    8  11  16  21
    2  14  23  27
    1   4   8   9
    8  9   12  17
    3  6   22  23
    7  14  22  27
    4   4   7   9
    1   6  18  19
    3  14  18  23
   10  10  23  27
    2   6   9  11
    6   6  17  19
    6  13  18  23
    3  16  24  29
    6   6   7  11
    6  10  15  19
    9  12  20  25
   11  12  24  29
    3   4  12  13
    4   5  20  21
   12  15  16  25
   12  16  21  29
    2   5  14  15
    4   8  19  21
    2   7  26  27 ];
# ------------------------------------------------------------------------------
# ---------------------------------------------- matrices and vectors of symbols
# ------------------------------------------------------------------------------
"""
    symbol_vector(s, indices) -> Vector{Symbol}

Construct a vector of symbolic names using the base string `s` and the supplied
indices.

For example, `symbol_vector("x", 1:3)` returns `[:x_1, :x_2, :x_3]`.
"""
function symbol_vector( s, indices )
    sstr = string(s)
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
    sstr = string(s)
    rows = collect(row_indices)
    cols = collect(col_indices)
    matrix = [Symbol("$(sstr)_{$(i),$(j)}") for i in rows, j in cols]
    return matrix
end
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
# ------------------------------------------------------------------------------
# -------------------------------------------------------------- GE, GJ problems
# ------------------------------------------------------------------------------
"""
    gen_gj_matrix(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false) -> pivot_cols, A

Generate a matrix suitable for Gauss-Jordan style exercises.
The returned matrix has rank `r`, controlled pivot placement, and exact integer
entries. `pivot_cols` records the pivot columns of the underlying echelon model.
"""
function gen_gj_matrix(m,n,r; maxint=3, pivot_in_first_col=true, has_zeros=false )
    M,pivot_cols=rref_matrix(m,n,r,maxint=maxint,pivot_in_first_col=pivot_in_first_col, has_zeros=has_zeros )

    s = ones( Int, n )
    s[pivot_cols] = rand( [-maxint:-1;1:maxint], r )

    A = unit_lower(m,maxint=maxint) * unit_lower(m,maxint=maxint)' * M * Diagonal(s)
    pivot_cols, A
end
# ------------------------------------------------------------------------------
"""
    gen_rhs(A, pivot_cols; maxint=3, num_rhs=1, has_zeros=false) -> X, B

Generate a random RHS matrix `B = A*X` consistent with given pivot columns.
"""
function gen_rhs( A, pivot_cols; maxint=3,num_rhs=1,has_zeros=false)
    rng = _int_range(maxint,has_zeros)
    X   = zeros(Int64, (size(A,2),num_rhs))

    X[pivot_cols,:] = rand( rng, (length(pivot_cols),num_rhs) )
    B = A*X
    X,B
end
# ------------------------------------------------------------------------------
# given the pivot locations, generate a particular solution of N integer entries, free variables set to zero
"""
    gen_particular_solution(pivot_cols, n; maxint=3, num_rhs=1) -> Matrix{Int}

Generate one or more particular solution vectors with prescribed pivot entries.
Only the pivot coordinates listed in `pivot_cols` are assigned nonzero random
integers; all free-variable coordinates are set to zero.
"""
function gen_particular_solution( pivot_cols, n; maxint=3, num_rhs=1 )
    X               = zeros(Int64, (n,num_rhs))
    X[pivot_cols,:] = rand( [-maxint:-1; 1:maxint], (length(pivot_cols),num_rhs) )
    X
end
# ------------------------------------------------------------------------------
"""
    gen_gj_pb(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false, num_rhs=1) -> A, X, B

Generate a consistent linear system `A * X = B` of rank `r`.
The right-hand side is produced from an exact integer solution matrix `X`.
"""
function gen_gj_pb(m,n,r;
        maxint=3, pivot_in_first_col=true, has_zeros=false, num_rhs=1 )
    _validate_rank_request("gen_gj_pb", m, n, r)
    pivot_cols,A = gen_gj_matrix(m,n,r;
                                 maxint=maxint, pivot_in_first_col=pivot_in_first_col, has_zeros=has_zeros )
    X,B=gen_rhs(A, pivot_cols; maxint=maxint,num_rhs=num_rhs,has_zeros=has_zeros)
    A,X,B
end
# ------------------------------------------------------------------------------
"""
    gen_gj_pb(m, n; maxint=3) -> A, X, B

Convenience method for the full-rank case `r = min(m, n)`.
"""
function gen_gj_pb(m,n; maxint=3)
    gen_gj_pb( m,n,min(m,n); maxint=maxint )
end
# ------------------------------------------------------------------------------
"""
    gen_inconsistent_gj_pb(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false, num_rhs=1) -> A, B

Generate an inconsistent linear system `A * x = B` of rank `r`.

The returned right-hand side columns are constructed outside the column space of
`A`, so each RHS column is inconsistent. This requires `r < m`, since a full
row-rank `m × n` matrix cannot produce an inconsistent system.
"""
function gen_inconsistent_gj_pb(m,n,r;
        maxint=3, pivot_in_first_col=true, has_zeros=false, num_rhs=1 )
    _validate_rank_request("gen_inconsistent_gj_pb", m, n, r)
    r < m || throw(ArgumentError("gen_inconsistent_gj_pb requires r < m so the system can be inconsistent"))

    M,pivot_cols = rref_matrix(m,n,r;
                               maxint=maxint,
                               pivot_in_first_col=pivot_in_first_col,
                               has_zeros=has_zeros )

    s = ones(Int, n)
    s[pivot_cols] = rand([-maxint:-1; 1:maxint], r)

    E = unit_lower(m,maxint=maxint) * unit_lower(m,maxint=maxint)'
    A = E * M * Diagonal(s)

    rng = _int_range(maxint, has_zeros)
    Y = zeros(Int64, m, num_rhs)
    if r > 0
        Y[1:r, :] = rand(rng, (r, num_rhs))
    end
    # A column of E*Y is inconsistent exactly when Y has a nonzero entry below
    # the rank rows of M, since the corresponding rows of M are zero.
    Y[r+1:m, :] = rand([-maxint:-1; 1:maxint], (m-r, num_rhs))
    B = E * Y

    A, B
end
# ------------------------------------------------------------------------------
"""
    invert_unit_lower(L) -> Matrix

Return the inverse of a unit lower-triangular matrix `L`.
The result has the same element type as `L`.
"""
function invert_unit_lower(L)
    n = size(L, 1)
    L_inv = Matrix{eltype(L)}(I, n, n)

    for j in n-1:-1:1 # current column of L_inv to update
        for k in j+1:n   # Use these columns of L_inv
            for i in k:n # each affected row
                L_inv[i, j] -= L[k, j] * L_inv[i, k]
            end
        end
    end
    return L_inv
end
# ------------------------------------------------------------------------------
"""
    gen_inv_pb(n; maxint=3) -> A, A_inv

Generate an invertible `n × n` integer matrix together with its exact inverse.
The construction is unimodular, so `A_inv` also has integer entries.
"""
function gen_inv_pb(n; maxint=3)
    # create an invertible matix problem of size n x n
    # with maxint=2, this works for n <= 15 or so
    e1 = unit_lower( n,n, maxint=maxint )
    e2 = unit_lower( n,n, maxint=maxint )
    A  = e1*e2'

    # A is unimodular: unit-lower factors and their transpose have determinant 1,
    # so the inverse remains integral.
    #A_inv = invert_unit_lower(e2)'*invert_unit_lower(e1)
    A_inv = Int.(inv(Rational{Int}.(A)))
    A, A_inv
end
# ------------------------------------------------------------------------------
"""
    gen_ldlt_pb(m; maxint=3, rank=nothing, squares=false) -> L, D, A

Generate an exact symmetric matrix in `L * D * L'` form.
When `rank` is provided, trailing diagonal entries of `D` are set to zero to
control the rank. When `squares=true`, the nonzero diagonal entries of `D` are
chosen from perfect squares.
"""
function gen_ldlt_pb(m;maxint=3,rank=nothing, squares = false)
    L   = unit_lower(m,maxint=maxint) 
    p   =  squares ? (1:maxint).^2 : 1:maxint
    if rank !== nothing
        0 <= rank <= m || throw(ArgumentError("rank must satisfy 0 <= rank <= m"))
        pivots = [rand( p, rank); zeros(Int, m-rank)]
        D   = Diagonal( pivots )
    else
        D   = Diagonal( rand( p, m))
    end

    A = L * D * L'
    L, D, A
end
# ------------------------------------------------------------------------------
"""
    gen_lu_pb(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false) -> pivot_cols, L, U, A

Generate an exact LU factorization exercise with `A = L * U`.
`U` is an `m × n` row-echelon matrix of rank `r`, `L` is unit lower triangular,
and `pivot_cols` records the pivot columns of `U`.
"""
function gen_lu_pb(m,n,r;maxint=3,pivot_in_first_col=true, has_zeros=false)
    _validate_rank_request("gen_lu_pb", m, n, r)
    U,pivot_cols = ref_matrix(m,n,r,maxint=maxint,pivot_in_first_col=pivot_in_first_col, has_zeros=has_zeros )
    L   = unit_lower(m,maxint=maxint)

    A = L * U
    pivot_cols, L, U, A
end
# ------------------------------------------------------------------------------
"""
    _plu_dependent_positions(m, r; nswaps=nothing) -> Vector{Int}

Choose the sorted insertion positions used by `gen_plu_pb` for dependent rows
that will be moved upward among the first `r` pivot rows.
"""
function _denominator_lcm(x)
    if x isa Rational
        return denominator(x)
    elseif x isa Complex{<:Rational}
        return lcm(denominator(real(x)), denominator(imag(x)))
    else
        return 1
    end
end

function _factor_out_denominator(A::AbstractArray)
    isempty(A) && return 1, copy(A)
    d = foldl(lcm, (_denominator_lcm(x) for x in A); init=1)
    return d, d .* A
end

function _plu_dependent_positions(m, r; nswaps=nothing)
    maxswaps = min(m - r, max(r - 1, 0))
    if nswaps === nothing
        k = maxswaps > 0 ? 1 : 0
    else
        nswaps >= 0 || throw(ArgumentError("nswaps must be nonnegative"))
        k = min(nswaps, maxswaps)
    end
    k == 0 && return Int[]

    positions = shuffle!(collect(2:r+k-1))[1:k]
    sort!(positions)
end

"""
    _plu_base_lower_factor(m, r, dependent_positions; maxint=3) -> Matrix{Int}

Build the dense base lower-triangular factor `L0` used by `gen_plu_pb`.
For each planned inserted dependent row, the entries in columns `1:k` are
resampled and columns `k+1:r` are zeroed so that the row depends only on the
first `k` pivot rows, while the lower-triangular tail to the right of the rank
block remains random.
"""
function _plu_base_lower_factor(m, r, dependent_positions; maxint=3)
    L = unit_lower(m; maxint=maxint)
    for (j, pos) in enumerate(dependent_positions)
        dep_row = r + j
        npivot_above = pos - j
        npivot_above > 0 || continue
        coeffs = rand([-maxint:-1; 1:maxint], npivot_above)
        L[dep_row, 1:npivot_above] .= coeffs
        if npivot_above < r
            L[dep_row, npivot_above+1:r] .= 0
        end
    end
    L
end

"""
    _plu_row_order(m, r, dependent_positions) -> Vector{Int}

Return the row order used to form the final matrix `A = P * L * U`.
Each dependent row `r + j` is inserted at the requested position in
`dependent_positions`, and the remaining pivot/unused rows keep their relative
order.
"""
function _plu_row_order(m, r, dependent_positions)
    k = length(dependent_positions)
    row_order = Vector{Int}(undef, m)
    dep_lookup = Set(dependent_positions)
    pivot_row = 1
    dep_row = r + 1

    for pos in 1:(r+k)
        if pos in dep_lookup
            row_order[pos] = dep_row
            dep_row += 1
        else
            row_order[pos] = pivot_row
            pivot_row += 1
        end
    end

    for pos in r+k+1:m
        row_order[pos] = dep_row
        dep_row += 1
    end

    row_order
end

"""
    _gen_plu_from_factors(U, r; maxint=3, nswaps=nothing, return_schedule=false) -> P, L, A
    _gen_plu_from_factors(U, r; maxint=3, nswaps=nothing, return_schedule=true) -> P, L, A, dependent_positions

Internal helper for `gen_plu_pb`. Keep the echelon factor `U` fixed, start from
a dense base unit lower-triangular factor `L0`, zero only the minimum entries
needed in the rank block so selected lower rows depend on controlled prefixes of
the pivot rows, then build a permutation `P` that inserts those dependent rows
upward. The final matrix is `A = P * L * U`.
"""
function _gen_plu_from_factors(U, r; maxint=3, nswaps=nothing, return_schedule=false)
    m = size(U, 1)
    dependent_positions = _plu_dependent_positions(m, r; nswaps=nswaps)
    L = _plu_base_lower_factor(m, r, dependent_positions; maxint=maxint)
    row_order = _plu_row_order(m, r, dependent_positions)
    P = [j == row_order[i] ? 1 : 0 for i in 1:m, j in 1:m]
    A = P * L * U

    if return_schedule
        return P, L, A, dependent_positions
    end
    return P, L, A
end

 # ------------------------------------------------------------------------------
"""
    gen_plu_pb(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false, nswaps=nothing) -> pivot_cols, P, L, U, A

Generate an exact PLU factorization exercise with `A = P * L * U`.
`U` is the canonical row-echelon factor. `L` first creates dependent rows from
the zero rows below rank `r`, and `P` then interleaves those dependent rows
among the pivot rows so standard elimination encounters forced row exchanges.
When `nswaps` is provided, the construction inserts up to
`min(nswaps, m-r, r-1)` such dependent rows; when omitted, it uses one swap
whenever that is possible. In particular, when `r == m` there are no zero rows
available for this construction, so no forced dependent-row insertions occur
and the returned permutation may be the identity.
"""
function gen_plu_pb(m,n,r;maxint=3,pivot_in_first_col=true, has_zeros=false, nswaps=nothing)
    _validate_rank_request("gen_plu_pb", m, n, r)
    U, pivot_cols = ref_matrix(m,n,r; maxint=maxint, pivot_in_first_col=pivot_in_first_col, has_zeros=has_zeros)
    P, L, A = _gen_plu_from_factors(U, r; maxint=maxint, nswaps=nswaps)
    pivot_cols, P, L, U, A
end
# ------------------------------------------------------------------------------
# ---------------------------------------------------------- orthogonal matrices
# ------------------------------------------------------------------------------
"""
    W_2_matrix() -> d, W

Return a `2 × 2` exact matrix `W` built from a Pythagorean triple such that
`W' * W = d^2 * I`.
"""
function W_2_matrix()
    a,b,c = PythagoreanNumberTriplets[ rand(1:size(PythagoreanNumberTriplets,1)), : ]
    c,[ a -b; b a]
end
# ------------------------------------------------------------------------------
"""
    Q_2_matrix() -> Matrix{Rational{Int}}

Return the exact orthogonal matrix obtained by normalizing `W_2_matrix()`.
"""
function Q_2_matrix()
    c,W = W_2_matrix()
    W // c
end
# ------------------------------------------------------------------------------
"""
    W_3_matrix(; maxint=3) -> d, W

Return a structured `3 × 3` exact matrix whose Gram matrix is diagonal.
The leading `2 × 2` block comes from a Pythagorean triple and the remaining
entry is an independent small integer.
"""
function W_3_matrix(; maxint=3)
    a,b,c = PythagoreanNumberTriplets[ rand(1:size(PythagoreanNumberTriplets,1)), : ]
    # The Pythagorean scale c governs the 2×2 rotation block only; the third
    # row/column uses an independent small integer so A' * A remains diagonal.
    A = [ a -b 0
          b  a 0
          0  0 rand( [-maxint:-1; 1:maxint]) ]
    A = A[shuffle(1:3),:]
    c,A[ :, shuffle(1:3)]
end
# ------------------------------------------------------------------------------
"""
    Q_3_matrix() -> Matrix{Rational{Int}}

Return the exact orthogonal matrix obtained by normalizing the `W_3_matrix`
construction.
"""
function Q_3_matrix()
    a,b,c = PythagoreanNumberTriplets[ rand(1:size(PythagoreanNumberTriplets,1)), : ]
    A = [ a//c -b//c  0
          b//c  a//c  0
             0     0  1 ]
    A = A[shuffle(1:3),:]
    A[ :, shuffle(1:3)]
end
# ------------------------------------------------------------------------------
# the following matrix has a block structure
"""
    Q_4_blocks() -> Matrix{Rational{Int}}

Return a `4 × 4` exact orthogonal matrix formed from two independent
Pythagorean `2 × 2` rotation blocks, followed by row and column shuffling.
"""
function Q_4_blocks()
    a1,b1,c1 = PythagoreanNumberTriplets[ rand(1:size(PythagoreanNumberTriplets,1)), : ]
    a2,b2,c2 = PythagoreanNumberTriplets[ rand(1:size(PythagoreanNumberTriplets,1)), : ]

    A = [ a1//c1 -b1//c1             0            0
          b1//c1  a1//c1             0            0
          0                          0   a2//c2  b2//c2 
          0                          0  -b2//c2  a2//c2 ]

    A = A[shuffle(1:4), :]
    A[ :, shuffle(1:4)]
end
# ------------------------------------------------------------------------------
"""
    W_4_matrix() -> d, W

Return a structured `4 × 4` exact matrix with diagonal Gram matrix.
The normalization factor `d` satisfies that `(W // d)` is orthogonal.
"""
function W_4_matrix()
    a,b,c,d = PythagoreanNumberQuadruplets[ rand(1:size(PythagoreanNumberQuadruplets,1)), : ]
    p  = (a*a + b*b) * d*d
    a2 = -a*c* p
    a3 =  a*b* p
    a4 =  a*a* p

    den = gcd(gcd( a2, a3), a4 )
    a2 = Int( a2 // den)
    a3 = Int( a3 // den)
    a4 = Int( a4 // den)

    A = [ a -b -c   0
          b  a  0  a2
          c  0  a  a3 
          0  c -b  a4 ]
    A = A[shuffle(1:4), :]
    d, A[:, shuffle(1:4)]
end
# ------------------------------------------------------------------------------
"""
    Q_4_matrix() -> Matrix{Rational{Int}}

Return the exact orthogonal matrix obtained by normalizing `W_4_matrix()`.
"""
function Q_4_matrix()
    d,W = W_4_matrix()
    W//d
end
# ------------------------------------------------------------------------------
"""
    W_matrix(n; general=false) -> d, W

Return a structured exact matrix `W` together with a common denominator `d`
such that `(W // d)` is orthogonal.
For `n == 2`, `3`, or `4` and `general=false`, specialized families are used;
otherwise the result is derived from `Q_matrix`.
"""
function W_matrix(n; general=false)
  if general == false
    if     n == 2 return W_2_matrix()
    elseif n == 3 return W_3_matrix()
    elseif n == 4 return W_4_matrix()
    end
  end
  A = Q_matrix(n; general=general)
  _factor_out_denominator(A)
end
# ------------------------------------------------------------------------------
"""
    Q_matrix(n; maxint=3, with_zeros=false, general=false) -> Matrix

Return an exact orthogonal matrix.
For sizes `2`, `3`, and `4` with `general=false`, this dispatches to the
specialized `Q_k_matrix` families. Otherwise it applies the Cayley transform to
a random skew-symmetric matrix.
"""
function Q_matrix(n; maxint=3, with_zeros=false, general=false )
  if general == false
    if     n == 2 return Q_2_matrix()
    elseif n == 3 return Q_3_matrix()
    elseif n == 4 return Q_4_matrix()
    end
  end
  S=skew_symmetric_matrix(n,maxint=maxint, with_zeros=with_zeros)
  inv(S-(1//1)I(size(S,1))) * (S+1I(size(S,1)))
end
# ------------------------------------------------------------------------------
"""
    sparse_Q_matrix(n; maxint=3, with_zeros=false) -> Matrix

Construct an exact rational orthogonal matrix from block sizes `n`.
For each block size `m` in `n`, the function builds a skew-symmetric
matrix and applies the Cayley transform `(S - I)^(-1) * (S + I)` to get
an orthogonal block. Those blocks are assembled into a block-diagonal
matrix and then randomly permuted by rows and columns.

When `n` is a single integer, this is the one-block case.
For very small block sizes such as `(2,2)`, the number of distinct
outputs can be limited when `maxint` is small; increase `maxint` if
more variation is desired.
"""
function sparse_Q_matrix(n; maxint=3, with_zeros=false )
    sz = sum(n)
    A  = zeros(Rational{Int64},(sz,sz))
    i  = 1
    for m in n
        S = Rational{Int64}.( skew_symmetric_matrix(m; maxint=maxint, with_zeros=with_zeros ) )
        E = (1//1)I(m)
        F = inv( S - E ) * ( S + E )
        rng = i:i+m-1 |> collect
        A[rng,rng] = F
        i += m
    end

    A = A[shuffle(1:sz), :]
    A[ :, shuffle(1:sz)]
end
# ------------------------------------------------------------------------------
"""
    sparse_W_matrix(n) -> d, W

Return the denominator-factor representation of `sparse_Q_matrix(n)`.
The pair satisfies `(W // d) == sparse_Q_matrix(n)`.
"""
function sparse_W_matrix(n)
    A = sparse_Q_matrix(n)
    _factor_out_denominator(A)
end
# ------------------------------------------------------------------------------
# ---------------------------------------------------------------- Orthogonality
# ------------------------------------------------------------------------------
"""
    ca_projection_matrix(A) -> Matrix

Return the orthogonal projection matrix onto the column space of `A`,
computed as `A * inv(A' * A) * A'`.

This assumes that the columns of `A` are linearly independent so that
`A' * A` is invertible.
"""
function ca_projection_matrix(A)
    rank(A) == size(A, 2) || throw(ArgumentError("ca_projection_matrix requires linearly independent columns so that A' * A is invertible"))
    A*inv(A'A)*A'
end
# ------------------------------------------------------------------------------
"""
    gen_qr_problem(n; family=:auto, maxint=3) -> Matrix

Generate a QR exercise matrix from one of several orthogonal-seed families.

- `family=:hadamard` uses the existing Hadamard-based construction and requires
  an integer size supported by `Hadamard.hadamard`.
- `family=:pythagorean` uses the specialized 2×2, 3×3, or 4×4 `W_k_matrix` families.
- `family=:cayley` uses a dense Cayley-transform orthogonal seed.
- `family=:sparse` uses `sparse_Q_matrix`, treating `n` as block sizes.
- `family=:auto` chooses `:pythagorean` for `n == 2`, `3`, or `4`, otherwise
  tries `:hadamard` and falls back to `:cayley` for integer sizes or `:sparse`
  for block-size inputs.
"""
function gen_qr_problem(n; family=:auto, maxint=3)
    function total_size(n)
        n isa Integer && return n
        return sum(n)
    end

    if family == :pythagorean
        if n == 2
            _, W = W_2_matrix()
            return W * unit_lower(2, maxint=maxint)'
        elseif n == 3
            _, W = W_3_matrix(maxint=maxint)
            return W * unit_lower(3, maxint=maxint)'
        elseif n == 4
            _, W = W_4_matrix()
            return W * unit_lower(4, maxint=maxint)'
        end
        throw(ArgumentError("family=:pythagorean is only supported for n == 2, 3, or 4"))
    elseif family == :auto && n isa Integer && n in (2, 3, 4)
        return gen_qr_problem(n; family=:pythagorean, maxint=maxint)
    end

    Qseed = _orthogonal_matrix_family(n; family=family, maxint=maxint)
    return Qseed * lower(total_size(n), maxint=maxint)'
end
# ------------------------------------------------------------------------------
"""
    gen_qr_problem_3(; maxint=3, kwargs...) -> Matrix

Legacy convenience wrapper for `gen_qr_problem(3; ...)`.
"""
function gen_qr_problem_3(; maxint=3, kwargs...)
    return gen_qr_problem(3; maxint=maxint, kwargs...)
end
# ------------------------------------------------------------------------------
"""
    gen_qr_problem_4(; maxint=3, kwargs...) -> Matrix

Legacy convenience wrapper for `gen_qr_problem(4; ...)`.
"""
function gen_qr_problem_4(; maxint=3, kwargs...)
    return gen_qr_problem(4; maxint=maxint, kwargs...)
end
# ------------------------------------------------------------------------------
"""
    _orthogonal_matrix_family(n; family=:auto, maxint=3) -> Matrix

Internal helper that returns the exact orthogonal matrix family selected by the
`family` keyword for `gen_qr_problem` and `gen_svd_problem`.
"""
function _orthogonal_matrix_family(n; family=:auto, maxint=3)
    if family == :auto
        if n isa Integer && n in (2, 3, 4)
            return _orthogonal_matrix_family(n; family=:pythagorean, maxint=maxint)
        elseif n isa Integer
            try
                return _orthogonal_matrix_family(n; family=:hadamard, maxint=maxint)
            catch err
                if err isa ArgumentError
                    return _orthogonal_matrix_family(n; family=:cayley, maxint=maxint)
                end
                rethrow()
            end
        else
            return _orthogonal_matrix_family(n; family=:sparse, maxint=maxint)
        end
    elseif family == :pythagorean
        n isa Integer || throw(ArgumentError("family=:pythagorean requires an integer size n"))
        if n == 2
            return Q_2_matrix()
        elseif n == 3
            return Q_3_matrix()
        elseif n == 4
            return Q_4_matrix()
        end
        throw(ArgumentError("family=:pythagorean is only supported for n == 2, 3, or 4"))
    elseif family == :hadamard
        n isa Integer || throw(ArgumentError("family=:hadamard requires an integer size n"))
        H = _ensure_hadamard()
        had = try
            Base.invokelatest(H.hadamard, n)
        catch err
            if err isa ArgumentError
                throw(ArgumentError("family=:hadamard requires a size supported by Hadamard.hadamard"))
            end
            rethrow()
        end
        d = isqrt(n)
        d * d == n || throw(ArgumentError("family=:hadamard requires n to have an integer square root for exact orthogonal factors"))
        Rational{Int64}.(had) ./ d
    elseif family == :cayley
        n isa Integer || throw(ArgumentError("family=:cayley requires an integer size n"))
        Q_matrix(n; maxint=maxint, general=true)
    elseif family == :sparse
        sparse_Q_matrix(n; maxint=maxint)
    else
        throw(ArgumentError("unknown orthogonal matrix family: $family"))
    end
end
# ------------------------------------------------------------------------------
# ---------------------------------------------------------------- Eigenproblems
# ------------------------------------------------------------------------------
"""
    _eigenvalues_vector(e_vals) -> Vector

Normalize a vector-like eigenvalue input into a one-dimensional Julia vector.
Accepts vectors and `1 × n` / `n × 1` matrices.
"""
function _eigenvalues_vector(e_vals)
    if e_vals isa AbstractVector
        return collect(e_vals)
    elseif e_vals isa AbstractMatrix
        min(size(e_vals)...) == 1 || throw(ArgumentError("e_vals must be a vector or a 1×n/n×1 matrix"))
        return vec(e_vals)
    end
    throw(ArgumentError("e_vals must be a vector or a 1×n/n×1 matrix"))
end

"""
    gen_eigenproblem(e_vals; maxint=3) -> S, Λ, S_inv, A

Generate a diagonalizable exact eigenproblem with prescribed eigenvalues.
`e_vals` may be a vector, a `1 × n` matrix, or an `n × 1` matrix. The return
values satisfy `A = S * Λ * S_inv`.
"""
function gen_eigenproblem( e_vals; maxint=3 )
    vals = _eigenvalues_vector(e_vals)
    Λ = Diagonal(vals)
    S,S_inv = gen_inv_pb(length(vals), maxint=maxint )
    S,Λ,S_inv, S*Λ*S_inv
end
# ------------------------------------------------------------------------------
"""
    gen_cx_eigenproblem(evals_no_conj; maxint=1) -> S, Λ, S_inv, A

Generate a real matrix whose eigenstructure includes the supplied complex
eigenvalues. Each nonreal scalar is expanded to its real `2 × 2` companion
block, so the returned `Λ` is block diagonal and `A = S * Λ * S_inv`.
"""
function gen_cx_eigenproblem( evals_no_conj; maxint=1 )
    function construct_diagonal_blocks()
        t = typeof( real( evals_no_conj[1] ))
        function f(x)
            if imag(x) == zero( t )
                [x]
            else
                [real(x) -imag(x); imag(x) real(x)]
            end
        end
        blocks  = [f(x) for x in evals_no_conj ]
        sz = sum( (x->size(x,1)).(blocks))
        A = zeros( t, sz, sz)
        k = 1
        for b in blocks
            l = size(b,1)-1
            A[k:k+l, k:k+l] = b
            k += l+1
        end
        A
    end

    Λ       = construct_diagonal_blocks()
    S,S_inv = gen_inv_pb( size(Λ,1), maxint=maxint )
    S,Λ,S_inv, S*Λ*S_inv
end
# ------------------------------------------------------------------------------
"""
    gen_symmetric_eigenproblem(e_vals; maxint=5, with_zeros=false, general=true) -> Q, Λ, A

Generate a symmetric exact eigenproblem with prescribed eigenvalues.
The return values satisfy `A = Q * Λ * Q'`, where `Q` is orthogonal.
"""
function gen_symmetric_eigenproblem( e_vals; maxint=5, with_zeros=false, general=true )
    vals = _eigenvalues_vector(e_vals)
    S = Q_matrix( length(vals); maxint=maxint, with_zeros=with_zeros, general=general )
    Λ = Diagonal(vals)
    S, Λ, S * Λ * S'
end
# ------------------------------------------------------------------------------
"""
    gen_non_diagonalizable_eigenproblem(e_dup, e; maxint=4) -> Matrix

Generate a `3 × 3` exact matrix with a repeated eigenvalue `e_dup` and a
nontrivial Jordan block, so the matrix is not diagonalizable.
"""
function gen_non_diagonalizable_eigenproblem( e_dup, e; maxint=4 )
    # size 3x3 problem
    S,S_inv = gen_inv_pb(3, maxint=maxint )
    Λ = [e_dup 1 0; 0 e_dup 0; 0 0 e]
    S * Λ * S_inv
end
# ------------------------------------------------------------------------------
"""
    jordan_block(λ, k) -> Bidiagonal

Construct the `k × k` Jordan block with eigenvalue `λ`.
"""
function jordan_block(lambda,k)
    k > 0 || throw(ArgumentError("k must be positive"))
    J = Bidiagonal( fill(lambda,k), ones(typeof(lambda),k-1),:U)
end
# ------------------------------------------------------------------------------
"""
    jordan_form(j_blocks) -> Matrix

Assemble a block-diagonal Jordan matrix from a collection of Jordan blocks.
"""
function jordan_form( j_blocks )
    sz = sum([ size(b,1) for b in j_blocks ])
    A  = zeros( eltype( j_blocks[1]), sz, sz )
    i = 1
    for b in j_blocks
        sz_b = size(b,1)
        A[i:i+sz_b-1, i:i+sz_b-1] = b
        i += sz_b
    end
    A
end
# ------------------------------------------------------------------------------
"""
    gen_from_jordan_form(j_blocks; maxint=3) -> Matrix

Generate a matrix similar to the block-diagonal Jordan form built from
`j_blocks`, using an exact invertible change of basis.
"""
function gen_from_jordan_form( j_blocks; maxint=3 )
    A = jordan_form( j_blocks )
    S,S_inv = gen_inv_pb( size(A,1), maxint=maxint )
    S*A*S_inv
end
# ------------------------------------------------------------------------------
# Generate a degenerate matrix based on block sizes or (eigenvalue, size) pairs
"""
    gen_degenerate_matrix(block_descriptions...; maxint=3) -> P, J, P_inv, A

Generate a matrix similar to a Jordan matrix with repeated or nilpotent blocks.
Each block description is either an integer `n` for a nilpotent `n × n` block or
a tuple `(λ, n)` for a Jordan block of size `n` with eigenvalue `λ`.
The return values satisfy `A = P * J * P_inv`.
"""
function gen_degenerate_matrix(block_descriptions::Vararg{Any}; maxint=3)
    total_size = 0
    for desc in block_descriptions
        if isa(desc, Int)
            desc > 0 || throw(ArgumentError("Jordan block sizes must be positive"))
            total_size += desc  # Integer block size (nilpotent case)
        elseif isa(desc, Tuple) && length(desc) == 2
            desc[2] > 0 || throw(ArgumentError("Jordan block sizes must be positive"))
            total_size += desc[2]  # Tuple (block eigenvalue, size)
        else
            throw(ArgumentError("Each block description must be an integer or a tuple (λ, n)."))
        end
    end

    function block_value_type(desc)
        if desc isa Int
            return Int
        elseif desc isa Tuple && length(desc) == 2
            return typeof(desc[1])
        else
            throw(ArgumentError("Each block description must be an integer or a tuple (λ, n)."))
        end
    end

    T = foldl(promote_type, (block_value_type(desc) for desc in block_descriptions))
    J = zeros(T, total_size, total_size)
    current_row = 1

    for desc in block_descriptions
        if isa(desc, Int)                                 # Nilpotent Jordan block
            n = desc
            J[current_row:(current_row+n-1), current_row:(current_row+n-1)] .= jordan_block(0, n)
        elseif isa(desc, Tuple) && length(desc) == 2      # Degenerate Jordan block with eigenvalue
            λ,n = desc
            J[current_row:(current_row+n-1), current_row:(current_row+n-1)] .= jordan_block(λ, n)
        end
        current_row += desc isa Int ? desc : desc[2]
    end

    P, P_inv = gen_inv_pb(total_size, maxint=maxint)
    P, J, P_inv, P * J * P_inv
end
# ------------------------------------------------------------------------------
"""
    gen_svd_problem(m, n, σ; left_family=:sparse, right_family=:sparse, maxint=3) -> U, Σ, Vt, A

Generate an SVD problem with specified singular values.
Returns orthogonal factors `U`, `Vt`, the diagonal matrix `Σ`, and the product `A = U * Σ * Vt`.
The `left_family` and `right_family` keywords choose the orthogonal-matrix
construction used for the left and right factors independently.
"""
function gen_svd_problem(m,n,σ; left_family=:sparse, right_family=:sparse, maxint = 3)
    U  = _orthogonal_matrix_family(m; family=left_family, maxint=maxint)
    Vt = _orthogonal_matrix_family(n; family=right_family, maxint=maxint)
    m = sum(m); n=sum(n)
    Σ  = zeros(eltype(σ[1]), m,n)
    for i in 1:min( m, size(σ,1) )
        Σ[i,i] = σ[i]
    end
    U, Σ, Vt, U * Σ * Vt
end
# ==============================================================================
