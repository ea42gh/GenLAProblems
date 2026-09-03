# ------------------------------------------------------------------------------
"""
    invert_unit_lower(L) -> Matrix

Return the inverse of a unit lower-triangular matrix `L`.
The result has the same element type as `L`.
"""
function invert_unit_lower(L)
    L isa AbstractMatrix || throw(ArgumentError("invert_unit_lower requires a matrix"))
    size(L, 1) == size(L, 2) ||
        throw(ArgumentError("invert_unit_lower requires a square matrix"))
    n = size(L, 1)
    all(L[i, i] == 1 for i = 1:n) ||
        throw(ArgumentError("invert_unit_lower requires unit diagonal entries"))
    all(iszero(L[i, j]) for i = 1:n for j = (i+1):n) ||
        throw(ArgumentError("invert_unit_lower requires a lower-triangular matrix"))

    L_inv = Matrix{eltype(L)}(I, n, n)

    for j = (n-1):-1:1 # current column of L_inv to update
        for k = (j+1):n   # Use these columns of L_inv
            for i = k:n # each affected row
                L_inv[i, j] -= L[k, j] * L_inv[i, k]
            end
        end
    end
    return L_inv
end
# ------------------------------------------------------------------------------
"""
    gen_inv_pb(n; maxint=3, rng=Random.default_rng()) -> A, A_inv

Generate an invertible `n × n` integer matrix together with its exact inverse.
The construction is unimodular, so `A_inv` also has integer entries.
"""
function gen_inv_pb(n; maxint = 3, rng = Random.default_rng())
    _validate_dimension("gen_inv_pb", "n", n)
    # create an invertible matix problem of size n x n
    # with maxint=2, this works for n <= 15 or so
    e1 = unit_lower(n, n, maxint = maxint, rng = rng)
    e2 = unit_lower(n, n, maxint = maxint, rng = rng)
    A = e1 * e2'

    # A is unimodular: unit-lower factors and their transpose have determinant 1,
    # so the inverse remains integral.
    #A_inv = invert_unit_lower(e2)'*invert_unit_lower(e1)
    A_inv = Int.(inv(Rational{Int}.(A)))
    A, A_inv
end
# ------------------------------------------------------------------------------
"""
    gen_ldlt_pb(m; maxint=3, rank=nothing, squares=false, rng=Random.default_rng()) -> L, D, A

Generate an exact symmetric matrix in `L * D * L'` form.
When `rank` is provided, trailing diagonal entries of `D` are set to zero to
control the rank. When `squares=true`, the nonzero diagonal entries of `D` are
chosen from perfect squares.
"""
function gen_ldlt_pb(
    m;
    maxint = 3,
    rank = nothing,
    squares = false,
    rng = Random.default_rng(),
)
    _validate_dimension("gen_ldlt_pb", "m", m)
    _validate_maxint("gen_ldlt_pb", maxint)
    L = unit_lower(m, maxint = maxint, rng = rng)
    p = squares ? (1:maxint) .^ 2 : 1:maxint
    if rank !== nothing
        _validate_dimension("gen_ldlt_pb", "rank", rank)
        rank <= m || throw(ArgumentError("rank must satisfy 0 <= rank <= m"))
        rank == 0 ||
            maxint > 0 ||
            throw(ArgumentError("gen_ldlt_pb requires maxint > 0 when rank is positive"))
        pivots = rank == 0 ? zeros(Int, m) : [rand(rng, p, rank); zeros(Int, m - rank)]
        D = Diagonal(pivots)
    else
        maxint > 0 || throw(ArgumentError("gen_ldlt_pb requires maxint > 0 for full rank"))
        D = Diagonal(rand(rng, p, m))
    end

    A = L * D * L'
    L, D, A
end
# ------------------------------------------------------------------------------
"""
    gen_lu_pb(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false, rng=Random.default_rng()) -> pivot_cols, L, U, A

Generate an exact LU factorization exercise with `A = L * U`.
`U` is an `m × n` row-echelon matrix of rank `r`, `L` is unit lower triangular,
and `pivot_cols` records the pivot columns of `U`.
"""
function gen_lu_pb(
    m,
    n,
    r;
    maxint = 3,
    pivot_in_first_col = true,
    has_zeros = false,
    rng = Random.default_rng(),
)
    _validate_rank_request("gen_lu_pb", m, n, r)
    U, pivot_cols = ref_matrix(
        m,
        n,
        r,
        maxint = maxint,
        pivot_in_first_col = pivot_in_first_col,
        has_zeros = has_zeros,
        rng = rng,
    )
    L = unit_lower(m, maxint = maxint, rng = rng)

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
    d = foldl(lcm, (_denominator_lcm(x) for x in A); init = 1)
    return d, d .* A
end

function _plu_dependent_positions(m, r; nswaps = nothing, rng = Random.default_rng())
    maxswaps = min(m - r, max(r - 1, 0))
    if nswaps === nothing
        k = maxswaps > 0 ? 1 : 0
    else
        nswaps >= 0 || throw(ArgumentError("nswaps must be nonnegative"))
        k = min(nswaps, maxswaps)
    end
    k == 0 && return Int[]

    positions = shuffle(rng, collect(2:(r+k-1)))[1:k]
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
function _plu_base_lower_factor(
    m,
    r,
    dependent_positions;
    maxint = 3,
    rng = Random.default_rng(),
)
    L = unit_lower(m; maxint = maxint, rng = rng)
    for (j, pos) in enumerate(dependent_positions)
        dep_row = r + j
        npivot_above = pos - j
        npivot_above > 0 || continue
        coeffs = rand(rng, [(-maxint):-1; 1:maxint], npivot_above)
        L[dep_row, 1:npivot_above] .= coeffs
        if npivot_above < r
            L[dep_row, (npivot_above+1):r] .= 0
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

    for pos = 1:(r+k)
        if pos in dep_lookup
            row_order[pos] = dep_row
            dep_row += 1
        else
            row_order[pos] = pivot_row
            pivot_row += 1
        end
    end

    for pos = (r+k+1):m
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
function _gen_plu_from_factors(
    U,
    r;
    maxint = 3,
    nswaps = nothing,
    return_schedule = false,
    rng = Random.default_rng(),
)
    m = size(U, 1)
    dependent_positions = _plu_dependent_positions(m, r; nswaps = nswaps, rng = rng)
    L = _plu_base_lower_factor(m, r, dependent_positions; maxint = maxint, rng = rng)
    row_order = _plu_row_order(m, r, dependent_positions)
    P = [j == row_order[i] ? 1 : 0 for i = 1:m, j = 1:m]
    A = P * L * U

    if return_schedule
        return P, L, A, dependent_positions
    end
    return P, L, A
end

# ------------------------------------------------------------------------------
"""
    gen_plu_pb(m, n, r; maxint=3, pivot_in_first_col=true, has_zeros=false, nswaps=nothing, rng=Random.default_rng()) -> pivot_cols, P, L, U, A

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
function gen_plu_pb(
    m,
    n,
    r;
    maxint = 3,
    pivot_in_first_col = true,
    has_zeros = false,
    nswaps = nothing,
    rng = Random.default_rng(),
)
    _validate_rank_request("gen_plu_pb", m, n, r)
    U, pivot_cols = ref_matrix(
        m,
        n,
        r;
        maxint = maxint,
        pivot_in_first_col = pivot_in_first_col,
        has_zeros = has_zeros,
        rng = rng,
    )
    P, L, A = _gen_plu_from_factors(U, r; maxint = maxint, nswaps = nswaps, rng = rng)
    pivot_cols, P, L, U, A
end
