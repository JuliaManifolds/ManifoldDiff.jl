
@doc raw"""
    subgrad_distance(M, q, p[, c = 2; atol = 0])
    subgrad_distance!(M, X, q, p[, c = 2; atol = 0])

compute the subgradient of the distance (in place of `X`)

```math
f(p) = \frac{1}{c} d^c_{\mathcal M}(p, q)
```

to a fixed point `q` on the manifold `M` and `c` is an
integer. The subgradient reads

```math
\partial f(p) = -d_{\mathcal M}^{c-2}(p, q)\log_pq
```

for ``c\neq 1`` or ``p\neq  q``. Note that for the remaining case ``c=1``,
``p=q``, the function is not differentiable. In this case, the subgradient is given by a tangent vector at `p` with norm less than or equal to one.

# Optional

* `c` – (`2`) the exponent of the distance,  i.e. the default is the distance
* `atol` – (`0`) the tolerance to use when evaluating the distance between `p` and `q`.
"""
function subgrad_distance(M, q, p, c::Int = 2; atol = zero(eltype(first(p))))
    if c == 2
        return -log(M, p, q)
    elseif c == 1 && distance(M, q, p) ≤ atol
        return normal_cone_vector(M, p)
    else
        return -distance(M, p, q)^(c - 2) * log(M, p, q)
    end
end
function subgrad_distance!(M, X, q, p, c::Int = 2; atol = zero(eltype(first(p))))
    log!(M, X, p, q)
    if c == 2
        X .*= -one(eltype(first(X)))
    elseif c == 1 && distance(M, q, p) ≤ atol
        normal_cone_vector!(M, X, p)
    else
        X .*= -distance(M, p, q)^(c - 2)
    end
    return X
end
function normal_cone_vector(M, p)
    Y = rand(M; vector_at = p)
    if norm(M, p, Y) > one(eltype(first(Y)))
        Y ./= norm(M, p, Y)
        Y .*= rand()
    end
    return Y
end
function normal_cone_vector!(M, Y, p)
    ManifoldsBase.rand!(M, Y; vector_at = p)
    if norm(M, p, Y) > one(eltype(first(Y)))
        Y ./= norm(M, p, Y)
        Y .*= rand()
    end
    return Y
end
