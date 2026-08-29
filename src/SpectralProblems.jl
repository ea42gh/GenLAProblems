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
        min(size(e_vals)...) == 1 ||
            throw(ArgumentError("e_vals must be a vector or a 1×n/n×1 matrix"))
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
function gen_eigenproblem(e_vals; maxint = 3)
    vals = _eigenvalues_vector(e_vals)
    Λ = Diagonal(vals)
    S, S_inv = gen_inv_pb(length(vals), maxint = maxint)
    S, Λ, S_inv, S * Λ * S_inv
end
# ------------------------------------------------------------------------------
"""
    gen_cx_eigenproblem(evals_no_conj; maxint=1) -> S, Λ, S_inv, A

Generate a real matrix whose eigenstructure includes the supplied complex
eigenvalues. Each nonreal scalar is expanded to its real `2 × 2` companion
block, so the returned `Λ` is block diagonal and `A = S * Λ * S_inv`.
"""
function gen_cx_eigenproblem(evals_no_conj; maxint = 1)
    function construct_diagonal_blocks()
        t = typeof(real(evals_no_conj[1]))
        function f(x)
            if imag(x) == zero(t)
                [x]
            else
                [real(x) -imag(x); imag(x) real(x)]
            end
        end
        blocks = [f(x) for x in evals_no_conj]
        sz = sum((x -> size(x, 1)).(blocks))
        A = zeros(t, sz, sz)
        k = 1
        for b in blocks
            l = size(b, 1) - 1
            A[k:k+l, k:k+l] = b
            k += l + 1
        end
        A
    end

    Λ = construct_diagonal_blocks()
    S, S_inv = gen_inv_pb(size(Λ, 1), maxint = maxint)
    S, Λ, S_inv, S * Λ * S_inv
end
# ------------------------------------------------------------------------------
"""
    gen_symmetric_eigenproblem(e_vals; maxint=5, with_zeros=false, general=true) -> Q, Λ, A

Generate a symmetric exact eigenproblem with prescribed eigenvalues.
The return values satisfy `A = Q * Λ * Q'`, where `Q` is orthogonal.
"""
function gen_symmetric_eigenproblem(e_vals; maxint = 5, with_zeros = false, general = true)
    vals = _eigenvalues_vector(e_vals)
    S = Q_matrix(length(vals); maxint = maxint, with_zeros = with_zeros, general = general)
    Λ = Diagonal(vals)
    S, Λ, S * Λ * S'
end
# ------------------------------------------------------------------------------
"""
    gen_non_diagonalizable_eigenproblem(e_dup, e; maxint=4) -> Matrix

Generate a `3 × 3` exact matrix with a repeated eigenvalue `e_dup` and a
nontrivial Jordan block, so the matrix is not diagonalizable.
"""
function gen_non_diagonalizable_eigenproblem(e_dup, e; maxint = 4)
    # size 3x3 problem
    S, S_inv = gen_inv_pb(3, maxint = maxint)
    Λ = [e_dup 1 0; 0 e_dup 0; 0 0 e]
    S * Λ * S_inv
end
# ------------------------------------------------------------------------------
"""
    jordan_block(λ, k) -> Bidiagonal

Construct the `k × k` Jordan block with eigenvalue `λ`.
"""
function jordan_block(lambda, k)
    k > 0 || throw(ArgumentError("k must be positive"))
    J = Bidiagonal(fill(lambda, k), ones(typeof(lambda), k - 1), :U)
end
# ------------------------------------------------------------------------------
"""
    jordan_form(j_blocks) -> Matrix

Assemble a block-diagonal Jordan matrix from a collection of Jordan blocks.
"""
function jordan_form(j_blocks)
    sz = sum([size(b, 1) for b in j_blocks])
    A = zeros(eltype(j_blocks[1]), sz, sz)
    i = 1
    for b in j_blocks
        sz_b = size(b, 1)
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
function gen_from_jordan_form(j_blocks; maxint = 3)
    A = jordan_form(j_blocks)
    S, S_inv = gen_inv_pb(size(A, 1), maxint = maxint)
    S * A * S_inv
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
function gen_degenerate_matrix(block_descriptions::Vararg{Any}; maxint = 3)
    total_size = 0
    for desc in block_descriptions
        if isa(desc, Int)
            desc > 0 || throw(ArgumentError("Jordan block sizes must be positive"))
            total_size += desc  # Integer block size (nilpotent case)
        elseif isa(desc, Tuple) && length(desc) == 2
            desc[2] > 0 || throw(ArgumentError("Jordan block sizes must be positive"))
            total_size += desc[2]  # Tuple (block eigenvalue, size)
        else
            throw(
                ArgumentError(
                    "Each block description must be an integer or a tuple (λ, n).",
                ),
            )
        end
    end

    function block_value_type(desc)
        if desc isa Int
            return Int
        elseif desc isa Tuple && length(desc) == 2
            return typeof(desc[1])
        else
            throw(
                ArgumentError(
                    "Each block description must be an integer or a tuple (λ, n).",
                ),
            )
        end
    end

    T = foldl(promote_type, (block_value_type(desc) for desc in block_descriptions))
    J = zeros(T, total_size, total_size)
    current_row = 1

    for desc in block_descriptions
        if isa(desc, Int)                                 # Nilpotent Jordan block
            n = desc
            J[current_row:(current_row+n-1), current_row:(current_row+n-1)] .=
                jordan_block(0, n)
        elseif isa(desc, Tuple) && length(desc) == 2      # Degenerate Jordan block with eigenvalue
            λ, n = desc
            J[current_row:(current_row+n-1), current_row:(current_row+n-1)] .=
                jordan_block(λ, n)
        end
        current_row += desc isa Int ? desc : desc[2]
    end

    P, P_inv = gen_inv_pb(total_size, maxint = maxint)
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
function gen_svd_problem(
    m,
    n,
    σ;
    left_family = :sparse,
    right_family = :sparse,
    maxint = 3,
    rng = Random.default_rng(),
)
    U = _orthogonal_matrix_family(m; family = left_family, maxint = maxint, rng = rng)
    Vt = _orthogonal_matrix_family(n; family = right_family, maxint = maxint, rng = rng)
    m = sum(m)
    n = sum(n)
    Σ = zeros(eltype(σ[1]), m, n)
    for i = 1:min(m, size(σ, 1))
        Σ[i, i] = σ[i]
    end
    U, Σ, Vt, U * Σ * Vt
end
# ==============================================================================
