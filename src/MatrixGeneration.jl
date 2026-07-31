#using LinearAlgebra, Latexify, Random

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
