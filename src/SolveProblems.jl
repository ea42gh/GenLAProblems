# ------------------------------------------------------------------------------
# ---------------------------------------------------------------- QR algoorithm
# ------------------------------------------------------------------------------
"""
    gram_schmidt_w(A) -> Matrix

Naive Gram-Schmidt producing an integer matrix of orthogonalized columns.
"""
function gram_schmidt_w(A)
    W = Array{Rational{Int64}}(undef, size(A))
    N = size(A,2)
    for j=1:N
        v_j = Rational.(A[:,j])
        for k=1:j-1
            v_j = v_j - (dot(W[:,k], A[:,j]) / dot(W[:,k], W[:,k]) ) .* W[:,k]
        end
        lcm_den = reduce((x, y) -> lcm(x, denominator(y)), v_j, init=1)
        tmp = lcm_den .* v_j
        tmp_num = numerator.(tmp)
        d = reduce(gcd, tmp_num, init=tmp_num[1])
        W[:, j] = tmp ./ d
    end
    W
end
# ------------------------------------------------------------------------------
"""
    normalize_columns(int_W) -> Matrix

Normalize integer columns, introducing symbolic square roots when needed.
"""
function normalize_columns( int_W )
    norms_squared = [dot(view(int_W, :, i), view(int_W, :, i)) for i in 1:size(int_W, 2)]
    norms         = []
    for norm_squared in norms_squared
        if norm_squared isa Rational
            if denominator(norm_squared) == 1
                sz = isqrt(numerator(norm_squared))
                if sz^2 == numerator(norm_squared)
                    push!(norms, sz)
                else
                    push!(norms, _ensure_symbolics().sqrt(norm_squared))
                end
            else
                push!(norms, _ensure_symbolics().sqrt(norm_squared))
            end
        else
            sz = isqrt(norm_squared)
            if sz^2 == norm_squared
                push!(norms, sz)
            else
                push!(norms, _ensure_symbolics().sqrt(norm_squared))
            end
        end
    end
    if all(x -> typeof(x) <: Integer, norms)
        norms = Int.(norms)
    end

    if eltype(norms) <: Integer
        Q = int_W .// norms'
    else
        Q = similar(int_W, eltype(norms))
        for i in 1:size(int_W, 2)
            Q[:, i] = view(int_W, :, i) ./ norms[i]'
        end
    end
    return Q
end
# ------------------------------------------------------------------------------
"""
    qr_layout(A) -> Any

Return a LaTeX-ready layout of the QR construction steps.
"""
function qr_layout(A)
    W = gram_schmidt_w(A)

    WtW  = Diagonal(W'W)
    WtA  = W'A
    S    =  ((x-> Rational{Int64}(round(sqrt(x)))).(WtW))^(-1)

    Qt = S * W'
    R  = S * WtA

    matrices =  [ [ nothing,  nothing,     A,          W ],
                  [ nothing,       W',   WtA,        WtW ],
                  [       S,       Qt,     R,    nothing ] ]

    to_latex( matrices )
end
# ------------------------------------------------------------------------------
"""
    gram_schmidt_stable(A::AbstractMatrix{T}; reorthogonalize=false) where T<:Number

Compute a stable Gram-Schmidt QR factorization, optionally reorthogonalizing.
"""
function gram_schmidt_stable(A::AbstractMatrix{T}; reorthogonalize=false) where T<:Number
    if T <: AbstractFloat || T <: Complex{<:AbstractFloat}
        Awork = Matrix{T}(A)
        m, n = size(Awork)
        Q = zeros(T, m, n)
        R = zeros(T, n, n)
        E = zeros(T, n, n)

        for j = 1:n
            v = Awork[:, j]

            if reorthogonalize
                for i = 1:j-1
                    E[i, j] = dot(Q[:, i], v)
                    v      -= E[i, j] * Q[:, i]
                end
            end

            R[j, j] = norm(v)
            Q[:, j] = v / R[j, j]

            for i = j+1:n
                E[j, i] = dot(Q[:, j], Awork[:, i])
                Awork[:, i] -= E[j, i] * Q[:, j]
            end
        end

        return Q, R
    end

    Awork = Matrix{Any}(A)
    m, n = size(Awork)
    Q = Matrix{Any}(undef, m, n)
    R = Matrix{Any}(undef, n, n)
    E = Matrix{Any}(undef, n, n)

    dot_any(x, y) = sum(conj(x[i]) * y[i] for i in eachindex(x, y))

    for j = 1:n
        v = Vector{Any}(Awork[:, j])

        if reorthogonalize
            for i = 1:j-1
                E[i, j] = dot_any(view(Q, :, i), v)
                v = v .- E[i, j] .* view(Q, :, i)
            end
        end

        norm_sq = dot_any(v, v)
        R[j, j] = _ensure_symbolics().sqrt(norm_sq)
        Q[:, j] = v ./ R[j, j]

        for i = j+1:n
            E[j, i] = dot_any(view(Q, :, j), Vector{Any}(Awork[:, i]))
            Awork[:, i] = Awork[:, i] .- E[j, i] .* Q[:, j]
        end
    end

    return Q, R
end
