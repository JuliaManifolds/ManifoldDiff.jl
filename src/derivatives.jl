@doc raw"""
    Y = geodesic_derivative(M, p, X, t::Number; γt = geodesic(M, p, X, t))
    geodesic_derivative!(M, Y, p, X, t::Number; γt = geodesic(M, p, X, t))

Evaluate the derivative of the geodesic ``γ(t)`` with ``γ_{p,X}(0) = p`` and ``\dot γ_{p,X}(0) = X`` at ``t``.
The formula reads

```math
\dot γ(t) = \mathcal P_{γ(t) \gets p} X
```

where $\mathcal P$ denotes the parallel transport.
This computation can also be done in-place of ``Y``.

# Optional Parameters

* `γt` – (`geodesic(M, p, X, t)`) the point on the geodesic at ``t``.
  This way if the point was computed earlier it can be resued here.
"""
function geodesic_derivative(M, p, X, t::Number; q = geodesic(M, p, X, t))
    return parallel_transport_to(M, p, X, q)
end
function geodesic_derivative!(M, Y, p, X, t::Number; q = geodesic(M, p, X, t))
    parallel_transport_to!(M, Y, p, X, q)
    return Y
end

@doc raw"""
    Y = shortest_geodesic_derivative(M, p, X, t::Number; γt = shortest_geodesic(M, p, q, t))
    shortest_geodesic_derivative!(M, Y, p, X, t::Number; γt = shortest_geodesic(M, p, q, t))

Evaluate the derivative of the shortest geodesic ``γ(t)`` with ``γ_{p,q}(0) = p`` and ``\dot γ_{p,q}(1) = q``
at ``t``.
The formula reads

```math
\dot γ(t) = \mathcal P_{γ(t) \gets p} \log_pq
```

where $\mathcal P$ denotes the parallel transport.
This computation can also be done in-place of ``Y``.

# Optional Parameters

* `γt = geodesic(M, p, X, t)` the point on the geodesic at ``t``.
  This way if the point was computed earlier it can be resued here.
"""
function shortest_geodesic_derivative(
    M,
    p,
    q,
    t::Number;
    γt = shortest_geodesic(M, p, q, t),
)
    return parallel_transport_to(M, p, log(M, p, q), γt)
end
function shortest_geodesic_derivative!(
    M,
    Y,
    p,
    q,
    t::Number;
    γt = shortest_geodesic(M, p, q, t),
)
    log!(M, Y, p, q)
    parallel_transport_to!(M, Y, p, Y, γt)
    return Y
end
