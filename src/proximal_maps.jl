@doc raw"""
    y = prox_distance(M, λ, p_data, p [, r=2])
    prox_distance!(M, q, λ, p_data, p [, r=2])

Compute the proximal map ``\operatorname{prox}_{λf}`` with
parameter λ of ``f(p) = \frac{1}{r}d_{\mathcal M}^r(p_data,p)``.
For the in-place variant the computation is done in place of `q`.

# Input
* `M`      a manifold `M`
* `λ`      the prox parameter, a positive real number.
* `p_data` a point on `M`.
* `p`      the argument of the proximal map
* `r`      (`2`) exponent of the distance.

# Output
* `q` – the result of the proximal map of ``f``

For more details see [WeinmannDemaretStorath:2014](@cite)
"""
function prox_distance(M::AbstractManifold, λ::Real, p_data, p, r::Int = 2)
    d = distance(M, p_data, p)
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
    return shortest_geodesic(M, p, p_data, t)
end
function prox_distance!(M::AbstractManifold, q, λ, p_data, p, r::Int = 2)
    d = distance(M, p_data, p)
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
    return shortest_geodesic!(M, q, p, p_data, t)
end
