
@doc raw"""
    grad_distance(M, q, p[, c=2])
    grad_distance!(M, X, q, p[, c=2])

compute the (sub)gradient of the distance (default: squared), in place of `X`.

```math
f(p) = \frac{1}{c} d^c_{\mathcal M}(p, q)
```

to a fixed point `q` on the manifold `M` and `c` is an
integer. The (sub-)gradient reads

```math
\operatorname{grad}f(p) = -d_{\mathcal M}^{c-2}(p, q)\log_pq
```

for ``c\neq 1`` or ``p\neq  q``. Note that for the remaining case ``c=1``,
``p=q``, the function is not differentiable. In this case, the function returns the
corresponding zero tangent vector, since this is an element of the subdifferential.

# Optional

* `c` – (`2`) the exponent of the distance,  i.e. the default is the squared
  distance
"""
function grad_distance(M, q, p, c::Int = 2)
    if c == 2
        return -log(M, p, q)
    elseif c == 1 && p == q
        return zero_vector(M, p)
    else
        return -distance(M, p, q)^(c - 2) * log(M, p, q)
    end
end
function grad_distance!(M, X, q, p, c::Int = 2)
    log!(M, X, p, q)
    if c == 2
        X .*= -one(eltype(X))
    elseif c == 1 && p == q
        X = zero_vector(M, p)
    else
        X .*= -distance(M, p, q)^(c - 2)
    end
    return X
end
