
@doc raw"""
    subgrad_distance(M, q, p[, c = 1, atol=eps(eltype(p))])
    subgrad_distance!(M, X, q, p[, c = 1, atol=eps(eltype(p))])

compute the subgradient of the distance (in place of `X`)

```math
f(p) = \frac{1}{c} d^c_{\mathcal M}(p, q)
```

to a fixed point `q` on the manifold `M` and `c` is an
integer. The subgradient reads

```math
\operatorname{grad}f(p) = -d_{\mathcal M}^{c-2}(p, q)\log_pq
```

for ``c\neq 1`` or ``p\neq  q``. Note that for the remaining case ``c=1``,
``p=q``, the function is not differentiable. In this case, the subgradient is given by a tangent vector at `p` with norm less than or equal to one.

# Optional

* `c` – (`1`) the exponent of the distance,  i.e. the default is the distance
* `atol` – (`eps(eltype(p))`) the tolerance to use when evaluating the distance between `p` and `q`.
"""
function subgrad_distance(M, q, p; c::Int = 1, atol = eps(eltype(p)))
    if c == 2
        return -log(M, p, q)
    elseif c == 1 && distance(M, q, p) ≤ atol
        X = rand(M; vector_at = p)
        (norm(M, p, X) > 1.0) && (X ./= norm(M, p, X))
        return X
    else
        return -distance(M, p, q)^(c - 2) * log(M, p, q)
    end
end
function subgrad_distance!(M, X, q, p; c::Int = 1, atol = eps(eltype(p)))
    log!(M, X, p, q)
    if c == 2
        X .*= -one(eltype(X))
    elseif c == 1 && distance(M, q, p) ≤ atol
        ManifoldsBase.rand!(M, X; vector_at = p)
        (norm(M, p, X) > 1.0) && (X ./= norm(M, p, X))
    else
        X .*= -distance(M, p, q)^(c - 2)
    end
    return X
end
