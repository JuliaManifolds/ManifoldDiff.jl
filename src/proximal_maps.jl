@doc raw"""
    y = prox_distance(M, λ, f, p [, r=2])
    prox_distance!(M, q, λ, f, p [, r=2])

compute the proximal map ``\operatorname{prox}_{λφ}`` with
parameter λ of ``φ(p) = \frac{1}{r}d_{\mathcal M}^r(f,p)``.
For the in-place variant the computation is done in place of `q`.

# Input
* `M` – a manifold `M`
* `λ` – the prox parameter
* `f` – a point ``f ∈ \mathcal M`` (the data)
* `p` – the argument of the proximal map
* `r` – (`2`) exponent of the distance.

# Output
* `q` – the result of the proximal map of ``φ``

For more details see [WeinmannDemaretStorath:2014](@cite)
"""
function prox_distance(M::AbstractManifold, λ, f, p, r::Int = 2)
    d = distance(M, f, p)
    if r == 2
        t = λ / (1 + λ)
    elseif r == 1
        t = (λ < d) ? λ / d : 1.0
    else
        throw(
            ErrorException(
                "Proximal Map of distance(M,f,x) not implemented for an exponent $(r) (requires 1 or 2)",
            ),
        )
    end
    return exp(M, p, log(M, p, f), t)
end
function prox_distance!(M::AbstractManifold, q, λ, f, p, r::Int = 2)
    d = distance(M, f, p)
    if r == 2
        t = λ / (1 + λ)
    elseif r == 1
        t = (λ < d) ? λ / d : 1.0
    else
        throw(
            ErrorException(
                "Proximal Map of distance(M,f,x) not implemented for an exponent $(r) (requires 1 or 2)",
            ),
        )
    end
    return exp!(M, q, p, log(M, p, f), t)
end
