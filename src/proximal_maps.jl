@doc raw"""
    y = prox_distance(M,λ,f,x [, p=2])
    prox_distance!(M, y, λ, f, x [, p=2])

compute the proximal map ``\operatorname{prox}_{λ\varphi}`` with
parameter λ of ``φ(x) = \frac{1}{p}d_{\mathcal M}^p(f,x)``.
For the mutating variant the computation is done in place of `y`.

# Input
* `M` – a manifold `M`
* `λ` – the prox parameter
* `f` – a point ``f ∈ \mathcal M`` (the data)
* `x` – the argument of the proximal map

# Optional argument
* `p` – (`2`) exponent of the distance.

# Output
* `y` – the result of the proximal map of ``φ``
"""
function prox_distance(M::AbstractManifold, λ, f, x, p::Int = 2)
    d = distance(M, f, x)
    if p == 2
        t = λ / (1 + λ)
    elseif p == 1
        t = (λ < d) ? λ / d : 1.0
    else
        throw(
            ErrorException(
                "Proximal Map of distance(M,f,x) not implemented for p=$(p) (requires p=1 or 2)",
            ),
        )
    end
    return exp(M, x, log(M, x, f), t)
end
function prox_distance!(M::AbstractManifold, y, λ, f, x, p::Int = 2)
    d = distance(M, f, x)
    if p == 2
        t = λ / (1 + λ)
    elseif p == 1
        t = (λ < d) ? λ / d : 1.0
    else
        throw(
            ErrorException(
                "Proximal Map of distance(M,f,x) not implemented for p=$(p) (requires p=1 or 2)",
            ),
        )
    end
    return exp!(M, y, x, log(M, x, f), t)
end
