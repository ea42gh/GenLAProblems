# ------------------------------------------------------------------------------
# ------------------------------------- matrices for use in GE and GJ algorithms
# ------------------------------------------------------------------------------

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
