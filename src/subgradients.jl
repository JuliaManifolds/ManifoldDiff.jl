
@doc raw"""
    subgrad_distance(M, q, p[, atol=eps(eltype(p))])
    subgrad_distance!(M, X, q, p[, atol=eps(eltype(p))])

compute the subgradient of the distance (in place of `X`)

```math
f(p) = d_{\mathcal M}(p, q)
```

to a fixed point `q` on the manifold `M`.
The subgradient at `p = q`is given by a norm one tangent vector at `p`.

# Optional

* `atol` – (`eps(eltype(p))`) the tolerance to use when evaluating the distance between `p` and `q`.
"""
function subgrad_distance(M, q, p; atol=eps(eltype(p)))
    @assert distance(M, q, p) ≤ atol ["Since p and q are different, use grad_distance."]
    X = rand(M; vector_at = p)
    (norm(M, p, X) > 1.0) && (X ./= norm(M, p, X))
    return X
end
function subgrad_distance!(M, X, q, p; atol=eps(eltype(p)))
    c = check_vector(M, p, X; atol = atol)
    if c === nothing
        @assert distance(M, q, p) ≤ atol ["Since p and q are different, use grad_distance."]
        (norm(M, p, X) > 1.0) && (X ./= norm(M, p, X))
        return X
    else
        println(c)
    end
end
