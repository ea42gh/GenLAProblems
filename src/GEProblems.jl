# ------------------------------------------------------------------------------
# -------------------------------------------------------------- GE, GJ problems
# ------------------------------------------------------------------------------
"""
    gen_gj_matrix(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false) -> pivot_cols, A

Generate a matrix suitable for Gauss-Jordan style exercises.
The returned matrix has rank `r`, controlled pivot placement, and exact integer
entries. `pivot_cols` records the pivot columns of the underlying echelon model.
"""
function gen_gj_matrix(m, n, r; maxint = 3, pivot_in_first_col = true, has_zeros = false)
    M, pivot_cols = rref_matrix(
        m,
        n,
        r,
        maxint = maxint,
        pivot_in_first_col = pivot_in_first_col,
        has_zeros = has_zeros,
    )

    s = ones(Int, n)
    s[pivot_cols] = rand([-maxint:-1; 1:maxint], r)

    A = unit_lower(m, maxint = maxint) * unit_lower(m, maxint = maxint)' * M * Diagonal(s)
    pivot_cols, A
end
# ------------------------------------------------------------------------------
"""
    gen_rhs(A, pivot_cols; maxint=3, num_rhs=1, has_zeros=false) -> X, B

Generate a random RHS matrix `B = A*X` consistent with given pivot columns.
"""
function gen_rhs(A, pivot_cols; maxint = 3, num_rhs = 1, has_zeros = false)
    rng = _int_range(maxint, has_zeros)
    X = zeros(Int64, (size(A, 2), num_rhs))

    X[pivot_cols, :] = rand(rng, (length(pivot_cols), num_rhs))
    B = A * X
    X, B
end
# ------------------------------------------------------------------------------
# given the pivot locations, generate a particular solution of N integer entries, free variables set to zero
"""
    gen_particular_solution(pivot_cols, n; maxint=3, num_rhs=1) -> Matrix{Int}

Generate one or more particular solution vectors with prescribed pivot entries.
Only the pivot coordinates listed in `pivot_cols` are assigned nonzero random
integers; all free-variable coordinates are set to zero.
"""
function gen_particular_solution(pivot_cols, n; maxint = 3, num_rhs = 1)
    X = zeros(Int64, (n, num_rhs))
    X[pivot_cols, :] = rand([-maxint:-1; 1:maxint], (length(pivot_cols), num_rhs))
    X
end
# ------------------------------------------------------------------------------
"""
    gen_gj_pb(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false, num_rhs=1) -> A, X, B

Generate a consistent linear system `A * X = B` of rank `r`.
The right-hand side is produced from an exact integer solution matrix `X`.
"""
function gen_gj_pb(
    m,
    n,
    r;
    maxint = 3,
    pivot_in_first_col = true,
    has_zeros = false,
    num_rhs = 1,
)
    _validate_rank_request("gen_gj_pb", m, n, r)
    pivot_cols, A = gen_gj_matrix(
        m,
        n,
        r;
        maxint = maxint,
        pivot_in_first_col = pivot_in_first_col,
        has_zeros = has_zeros,
    )
    X, B = gen_rhs(A, pivot_cols; maxint = maxint, num_rhs = num_rhs, has_zeros = has_zeros)
    A, X, B
end
# ------------------------------------------------------------------------------
"""
    gen_gj_pb(m, n; maxint=3) -> A, X, B

Convenience method for the full-rank case `r = min(m, n)`.
"""
function gen_gj_pb(m, n; maxint = 3)
    gen_gj_pb(m, n, min(m, n); maxint = maxint)
end
# ------------------------------------------------------------------------------
"""
    gen_inconsistent_gj_pb(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false, num_rhs=1) -> A, B

Generate an inconsistent linear system `A * x = B` of rank `r`.

The returned right-hand side columns are constructed outside the column space of
`A`, so each RHS column is inconsistent. This requires `r < m`, since a full
row-rank `m × n` matrix cannot produce an inconsistent system.
"""
function gen_inconsistent_gj_pb(
    m,
    n,
    r;
    maxint = 3,
    pivot_in_first_col = true,
    has_zeros = false,
    num_rhs = 1,
)
    _validate_rank_request("gen_inconsistent_gj_pb", m, n, r)
    r < m || throw(
        ArgumentError(
            "gen_inconsistent_gj_pb requires r < m so the system can be inconsistent",
        ),
    )

    M, pivot_cols = rref_matrix(
        m,
        n,
        r;
        maxint = maxint,
        pivot_in_first_col = pivot_in_first_col,
        has_zeros = has_zeros,
    )

    s = ones(Int, n)
    s[pivot_cols] = rand([-maxint:-1; 1:maxint], r)

    E = unit_lower(m, maxint = maxint) * unit_lower(m, maxint = maxint)'
    A = E * M * Diagonal(s)

    rng = _int_range(maxint, has_zeros)
    Y = zeros(Int64, m, num_rhs)
    if r > 0
        Y[1:r, :] = rand(rng, (r, num_rhs))
    end
    # A column of E*Y is inconsistent exactly when Y has a nonzero entry below
    # the rank rows of M, since the corresponding rows of M are zero.
    Y[r+1:m, :] = rand([-maxint:-1; 1:maxint], (m - r, num_rhs))
    B = E * Y

    A, B
end
