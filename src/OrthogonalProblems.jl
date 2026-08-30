# ------------------------------------------------------------------------------
# ---------------------------------------------------------- orthogonal matrices
# ------------------------------------------------------------------------------
"""
    W_2_matrix(; rng=Random.default_rng()) -> d, W

Return a `2 × 2` exact matrix `W` built from a Pythagorean triple such that
`W' * W = d^2 * I`.
"""
function W_2_matrix(; rng = Random.default_rng())
    a, b, c = PythagoreanNumberTriplets[rand(rng, 1:size(PythagoreanNumberTriplets, 1)), :]
    c, [a -b; b a]
end
# ------------------------------------------------------------------------------
"""
    Q_2_matrix(; rng=Random.default_rng()) -> Matrix{Rational{Int}}

Return the exact orthogonal matrix obtained by normalizing `W_2_matrix(rng = rng)`.
"""
function Q_2_matrix(; rng = Random.default_rng())
    c, W = W_2_matrix(rng = rng)
    W // c
end
# ------------------------------------------------------------------------------
"""
    W_3_matrix(; maxint=3, rng=Random.default_rng()) -> d, W

Return a structured `3 × 3` exact matrix whose Gram matrix is diagonal.
The leading `2 × 2` block comes from a Pythagorean triple and the remaining
entry is an independent small integer.
"""
function W_3_matrix(; maxint = 3, rng = Random.default_rng())
    _validate_positive_maxint("W_3_matrix", maxint)
    a, b, c = PythagoreanNumberTriplets[rand(rng, 1:size(PythagoreanNumberTriplets, 1)), :]
    # The Pythagorean scale c governs the 2×2 rotation block only; the third
    # row/column uses an independent small integer so A' * A remains diagonal.
    A = [
        a -b 0
        b a 0
        0 0 rand(rng, [(-maxint):-1; 1:maxint])
    ]
    A = A[shuffle(rng, 1:3), :]
    c, A[:, shuffle(rng, 1:3)]
end
# ------------------------------------------------------------------------------
"""
    Q_3_matrix(; rng=Random.default_rng()) -> Matrix{Rational{Int}}

Return the exact orthogonal matrix obtained by normalizing the `W_3_matrix`
construction.
"""
function Q_3_matrix(; rng = Random.default_rng())
    a, b, c = PythagoreanNumberTriplets[rand(rng, 1:size(PythagoreanNumberTriplets, 1)), :]
    A = [
        a//c -b//c 0
        b//c a//c 0
        0 0 1
    ]
    A = A[shuffle(rng, 1:3), :]
    A[:, shuffle(rng, 1:3)]
end
# ------------------------------------------------------------------------------
# the following matrix has a block structure
"""
    Q_4_blocks(; rng=Random.default_rng()) -> Matrix{Rational{Int}}

Return a `4 × 4` exact orthogonal matrix formed from two independent
Pythagorean `2 × 2` rotation blocks, followed by row and column shuffling.
"""
function Q_4_blocks(; rng = Random.default_rng())
    a1, b1, c1 =
        PythagoreanNumberTriplets[rand(rng, 1:size(PythagoreanNumberTriplets, 1)), :]
    a2, b2, c2 =
        PythagoreanNumberTriplets[rand(rng, 1:size(PythagoreanNumberTriplets, 1)), :]

    A = [
        a1//c1 -b1//c1 0 0
        b1//c1 a1//c1 0 0
        0 0 a2//c2 b2//c2
        0 0 -b2//c2 a2//c2
    ]

    A = A[shuffle(rng, 1:4), :]
    A[:, shuffle(rng, 1:4)]
end
# ------------------------------------------------------------------------------
"""
    W_4_matrix(; rng=Random.default_rng()) -> d, W

Return a structured `4 × 4` exact matrix with diagonal Gram matrix.
The normalization factor `d` satisfies that `(W // d)` is orthogonal.
"""
function W_4_matrix(; rng = Random.default_rng())
    a, b, c, d =
        PythagoreanNumberQuadruplets[rand(rng, 1:size(PythagoreanNumberQuadruplets, 1)), :]
    p = (a * a + b * b) * d * d
    a2 = -a * c * p
    a3 = a * b * p
    a4 = a * a * p

    den = gcd(gcd(a2, a3), a4)
    a2 = Int(a2 // den)
    a3 = Int(a3 // den)
    a4 = Int(a4 // den)

    A = [
        a -b -c 0
        b a 0 a2
        c 0 a a3
        0 c -b a4
    ]
    A = A[shuffle(rng, 1:4), :]
    d, A[:, shuffle(rng, 1:4)]
end
# ------------------------------------------------------------------------------
"""
    Q_4_matrix(; rng=Random.default_rng()) -> Matrix{Rational{Int}}

Return the exact orthogonal matrix obtained by normalizing `W_4_matrix(rng = rng)`.
"""
function Q_4_matrix(; rng = Random.default_rng())
    d, W = W_4_matrix(rng = rng)
    W // d
end
# ------------------------------------------------------------------------------
"""
    W_matrix(n; general=false, rng=Random.default_rng()) -> d, W

Return a structured exact matrix `W` together with a common denominator `d`
such that `(W // d)` is orthogonal.
For `n == 2`, `3`, or `4` and `general=false`, specialized families are used;
otherwise the result is derived from `Q_matrix`.
"""
function W_matrix(n; general = false, rng = Random.default_rng())
    if general == false
        if n == 2
            return W_2_matrix(rng = rng)
        elseif n == 3
            return W_3_matrix(rng = rng)
        elseif n == 4
            return W_4_matrix(rng = rng)
        end
    end
    A = Q_matrix(n; general = general, rng = rng)
    _factor_out_denominator(A)
end
# ------------------------------------------------------------------------------
"""
    Q_matrix(n; maxint=3, with_zeros=false, general=false, rng=Random.default_rng()) -> Matrix

Return an exact orthogonal matrix.
For sizes `2`, `3`, and `4` with `general=false`, this dispatches to the
specialized `Q_k_matrix` families. Otherwise it applies the Cayley transform to
a random skew-symmetric matrix.
"""
function Q_matrix(
    n;
    maxint = 3,
    with_zeros = false,
    general = false,
    rng = Random.default_rng(),
)
    _validate_maxint("Q_matrix", maxint)
    if general == false
        if n == 2
            return Q_2_matrix(rng = rng)
        elseif n == 3
            return Q_3_matrix(rng = rng)
        elseif n == 4
            return Q_4_matrix(rng = rng)
        end
    end
    S = skew_symmetric_matrix(n, maxint = maxint, with_zeros = with_zeros, rng = rng)
    inv(S - (1 // 1)I(size(S, 1))) * (S + 1I(size(S, 1)))
end
# ------------------------------------------------------------------------------
"""
    sparse_Q_matrix(n; maxint=3, with_zeros=false, rng=Random.default_rng()) -> Matrix

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
function sparse_Q_matrix(n; maxint = 3, with_zeros = false, rng = Random.default_rng())
    blocks = _validate_block_sizes("sparse_Q_matrix", n)
    sz = sum(blocks)
    A = zeros(Rational{Int64}, (sz, sz))
    i = 1
    for m in blocks
        S =
            Rational{
                Int64,
            }.(
                skew_symmetric_matrix(
                    m;
                    maxint = maxint,
                    with_zeros = with_zeros,
                    rng = rng,
                ),
            )
        E = (1 // 1)I(m)
        F = inv(S - E) * (S + E)
        inds = i:(i+m-1) |> collect
        A[inds, inds] = F
        i += m
    end

    A = A[shuffle(rng, 1:sz), :]
    A[:, shuffle(rng, 1:sz)]
end
# ------------------------------------------------------------------------------
"""
    sparse_W_matrix(n; rng=Random.default_rng()) -> d, W

Return the denominator-factor representation of `sparse_Q_matrix(n)`.
The pair satisfies `(W // d) == sparse_Q_matrix(n)`.
"""
function sparse_W_matrix(n; rng = Random.default_rng())
    A = sparse_Q_matrix(n; rng = rng)
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
    A isa AbstractMatrix ||
        throw(ArgumentError("ca_projection_matrix requires A to be a matrix"))
    rank(A) == size(A, 2) || throw(
        ArgumentError(
            "ca_projection_matrix requires linearly independent columns so that A' * A is invertible",
        ),
    )
    A * inv(A'A) * A'
end
# ------------------------------------------------------------------------------
"""
    gen_qr_problem(n; family=:auto, maxint=3, rng=Random.default_rng()) -> Matrix

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
function gen_qr_problem(n; family = :auto, maxint = 3, rng = Random.default_rng())
    n =
        n isa Integer ? (_validate_dimension("gen_qr_problem", "n", n; minimum = 1); n) :
        _validate_block_sizes("gen_qr_problem", n)

    function total_size(n)
        n isa Integer && return n
        return sum(n)
    end

    if family == :pythagorean
        if n == 2
            _, W = W_2_matrix(rng = rng)
            return W * unit_lower(2, maxint = maxint, rng = rng)'
        elseif n == 3
            _, W = W_3_matrix(maxint = maxint, rng = rng)
            return W * unit_lower(3, maxint = maxint, rng = rng)'
        elseif n == 4
            _, W = W_4_matrix(rng = rng)
            return W * unit_lower(4, maxint = maxint, rng = rng)'
        end
        throw(ArgumentError("family=:pythagorean is only supported for n == 2, 3, or 4"))
    elseif family == :auto && n isa Integer && n in (2, 3, 4)
        return gen_qr_problem(n; family = :pythagorean, maxint = maxint, rng = rng)
    end

    Qseed = _orthogonal_matrix_family(n; family = family, maxint = maxint, rng = rng)
    return Qseed * lower(total_size(n), maxint = maxint, rng = rng)'
end
"""
    _orthogonal_matrix_family(n; family=:auto, maxint=3) -> Matrix

Internal helper that returns the exact orthogonal matrix family selected by the
`family` keyword for `gen_qr_problem` and `gen_svd_problem`.
"""
function _orthogonal_matrix_family(
    n;
    family = :auto,
    maxint = 3,
    rng = Random.default_rng(),
)
    if family == :auto
        if n isa Integer && n in (2, 3, 4)
            return _orthogonal_matrix_family(
                n;
                family = :pythagorean,
                maxint = maxint,
                rng = rng,
            )
        elseif n isa Integer
            try
                return _orthogonal_matrix_family(
                    n;
                    family = :hadamard,
                    maxint = maxint,
                    rng = rng,
                )
            catch err
                if err isa ArgumentError
                    return _orthogonal_matrix_family(
                        n;
                        family = :cayley,
                        maxint = maxint,
                        rng = rng,
                    )
                end
                rethrow()
            end
        else
            return _orthogonal_matrix_family(
                n;
                family = :sparse,
                maxint = maxint,
                rng = rng,
            )
        end
    elseif family == :pythagorean
        n isa Integer ||
            throw(ArgumentError("family=:pythagorean requires an integer size n"))
        if n == 2
            return Q_2_matrix(rng = rng)
        elseif n == 3
            return Q_3_matrix(rng = rng)
        elseif n == 4
            return Q_4_matrix(rng = rng)
        end
        throw(ArgumentError("family=:pythagorean is only supported for n == 2, 3, or 4"))
    elseif family == :hadamard
        n isa Integer || throw(ArgumentError("family=:hadamard requires an integer size n"))
        had = try
            Base.invokelatest(Hadamard.hadamard, n)
        catch err
            if err isa ArgumentError
                throw(
                    ArgumentError(
                        "family=:hadamard requires a size supported by Hadamard.hadamard",
                    ),
                )
            end
            rethrow()
        end
        d = isqrt(n)
        d * d == n || throw(
            ArgumentError(
                "family=:hadamard requires n to have an integer square root for exact orthogonal factors",
            ),
        )
        Rational{Int64}.(had) ./ d
    elseif family == :cayley
        n isa Integer || throw(ArgumentError("family=:cayley requires an integer size n"))
        Q_matrix(n; maxint = maxint, general = true, rng = rng)
    elseif family == :sparse
        sparse_Q_matrix(n; maxint = maxint, rng = rng)
    else
        throw(ArgumentError("unknown orthogonal matrix family: $family"))
    end
end
