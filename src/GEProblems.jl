# ------------------------------------------------------------------------------
# -------------------------------------------------------------- GE, GJ problems
# ------------------------------------------------------------------------------
"""
    gen_gj_matrix(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false, rng=Random.default_rng()) -> pivot_cols, A

Generate a matrix suitable for Gauss-Jordan style exercises.
The returned matrix has rank `r`, controlled pivot placement, and exact integer
entries. `pivot_cols` records the pivot columns of the underlying echelon model.
"""
function gen_gj_matrix(
    m,
    n,
    r;
    maxint = 3,
    pivot_in_first_col = true,
    has_zeros = false,
    rng = Random.default_rng(),
)
    r == 0 || _validate_positive_maxint("gen_gj_matrix", maxint)

    M, pivot_cols = rref_matrix(
        m,
        n,
        r,
        maxint = maxint,
        pivot_in_first_col = pivot_in_first_col,
        has_zeros = has_zeros,
        rng = rng,
    )

    s = ones(Int, n)
    if r > 0
        s[pivot_cols] = rand(rng, [-maxint:-1; 1:maxint], r)
    end

    A =
        unit_lower(m, maxint = maxint, rng = rng) *
        unit_lower(m, maxint = maxint, rng = rng)' *
        M *
        Diagonal(s)
    pivot_cols, A
end
function _validate_num_rhs(fname, num_rhs)
    num_rhs isa Integer && num_rhs >= 0 ||
        throw(ArgumentError("$fname requires num_rhs to be a nonnegative integer"))
end

function _validate_pivot_columns(fname, pivot_cols, n)
    cols = pivot_cols isa Integer ? [pivot_cols] : collect(pivot_cols)
    all(x -> x isa Integer && 1 <= x <= n, cols) && length(unique(cols)) == length(cols) ||
        throw(ArgumentError("$fname requires unique pivot columns in 1:n"))
    cols
end
# ------------------------------------------------------------------------------
"""
    gen_rhs(A, pivot_cols; maxint=3, num_rhs=1, has_zeros=false, rng=Random.default_rng()) -> X, B

Generate a random RHS matrix `B = A*X` consistent with given pivot columns.
"""
function gen_rhs(
    A,
    pivot_cols;
    maxint = 3,
    num_rhs = 1,
    has_zeros = false,
    rng = Random.default_rng(),
)
    _validate_num_rhs("gen_rhs", num_rhs)
    _validate_maxint("gen_rhs", maxint)
    A isa AbstractMatrix || throw(ArgumentError("gen_rhs requires A to be a matrix"))
    pivot_cols = _validate_pivot_columns("gen_rhs", pivot_cols, size(A, 2))
    X = zeros(Int64, (size(A, 2), num_rhs))
    if !isempty(pivot_cols) && num_rhs > 0
        values = _int_range(maxint, has_zeros)
        X[pivot_cols, :] = rand(rng, values, (length(pivot_cols), num_rhs))
    end
    B = A * X
    X, B
end# ------------------------------------------------------------------------------
# given the pivot locations, generate a particular solution of N integer entries, free variables set to zero
"""
    gen_particular_solution(pivot_cols, n; maxint=3, num_rhs=1, rng=Random.default_rng()) -> Matrix{Int}

Generate one or more particular solution vectors with prescribed pivot entries.
Only the pivot coordinates listed in `pivot_cols` are assigned nonzero random
integers; all free-variable coordinates are set to zero.
"""
function gen_particular_solution(
    pivot_cols,
    n;
    maxint = 3,
    num_rhs = 1,
    rng = Random.default_rng(),
)
    _validate_maxint("gen_particular_solution", maxint)
    _validate_dimension("gen_particular_solution", "n", n)
    pivot_cols = _validate_pivot_columns("gen_particular_solution", pivot_cols, n)
    _validate_num_rhs("gen_particular_solution", num_rhs)
    X = zeros(Int64, (n, num_rhs))
    if !isempty(pivot_cols) && num_rhs > 0
        _validate_positive_maxint("gen_particular_solution", maxint)
        values = _int_range(maxint, false)
        X[pivot_cols, :] = rand(rng, values, (length(pivot_cols), num_rhs))
    end
    X
end# ------------------------------------------------------------------------------
"""
    gen_gj_pb(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false, num_rhs=1, rng=Random.default_rng()) -> A, X, B

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
    rng = Random.default_rng(),
)
    _validate_rank_request("gen_gj_pb", m, n, r)
    pivot_cols, A = gen_gj_matrix(
        m,
        n,
        r;
        maxint = maxint,
        pivot_in_first_col = pivot_in_first_col,
        has_zeros = has_zeros,
        rng = rng,
    )
    X, B = gen_rhs(
        A,
        pivot_cols;
        maxint = maxint,
        num_rhs = num_rhs,
        has_zeros = has_zeros,
        rng = rng,
    )
    A, X, B
end
# ------------------------------------------------------------------------------
"""
    gen_gj_pb(m, n; maxint=3, rng=Random.default_rng()) -> A, X, B

Convenience method for the full-rank case `r = min(m, n)`.
"""
function gen_gj_pb(m, n; maxint = 3, rng = Random.default_rng())
    gen_gj_pb(m, n, min(m, n); maxint = maxint, rng = rng)
end
# ------------------------------------------------------------------------------
"""
    gen_inconsistent_gj_pb(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false, num_rhs=1, rng=Random.default_rng()) -> A, B

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
    rng = Random.default_rng(),
)
    _validate_rank_request("gen_inconsistent_gj_pb", m, n, r)
    _validate_num_rhs("gen_inconsistent_gj_pb", num_rhs)
    (r == 0 && num_rhs == 0) || _validate_positive_maxint("gen_inconsistent_gj_pb", maxint)
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
        rng = rng,
    )

    s = ones(Int, n)
    if r > 0
        s[pivot_cols] = rand(rng, [-maxint:-1; 1:maxint], r)
    end

    E =
        unit_lower(m, maxint = maxint, rng = rng) *
        unit_lower(m, maxint = maxint, rng = rng)'
    A = E * M * Diagonal(s)

    Y = zeros(Int64, m, num_rhs)
    if num_rhs > 0
        values = _int_range(maxint, has_zeros)
        if r > 0
            Y[1:r, :] = rand(rng, values, (r, num_rhs))
        end
        # A column of E*Y is inconsistent exactly when Y has a nonzero entry below
        # the rank rows of M, since the corresponding rows of M are zero.
        Y[r+1:m, :] = rand(rng, [-maxint:-1; 1:maxint], (m - r, num_rhs))
    end
    B = E * Y

    A, B
end
