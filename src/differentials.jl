@doc raw"""
    Y = differential_shortest_geodesic_startpoint(M, p, q, t, X)
    differential_shortest_geodesic_startpoint!(M, Y, p, q, t, X)

Compute ``D_p γ(t;p,q)[η]`` (in place of `Y`).

# See also
[`differential_shortest_geodesic_endpoint`](@ref), [`jacobi_field`](@ref)
"""
function differential_shortest_geodesic_startpoint(M::AbstractManifold, p, q, t, X)
    return jacobi_field(M, p, q, t, X, βdifferential_shortest_geodesic_startpoint)
end
function differential_shortest_geodesic_startpoint!(M::AbstractManifold, Y, p, q, t, X)
    return jacobi_field!(M, Y, p, q, t, X, βdifferential_shortest_geodesic_startpoint)
end


@doc raw"""
    Y = differential_shortest_geodesic_endpoint(M, p, q, t, X)
    differential_shortest_geodesic_endpoint!(M, Y, p, q, t, X)

Compute ``D_qγ(t;p,q)[X]`` (in place of `Y`).

# See also
 [`differential_shortest_geodesic_startpoint`](@ref), [`jacobi_field`](@ref)
"""
function differential_shortest_geodesic_endpoint(M::AbstractManifold, p, q, t, X)
    return jacobi_field(M, q, p, 1 - t, X, βdifferential_shortest_geodesic_startpoint)
end
function differential_shortest_geodesic_endpoint!(M::AbstractManifold, Y, p, q, t, X)
    return jacobi_field!(M, Y, q, p, 1 - t, X, βdifferential_shortest_geodesic_startpoint)
end

@doc raw"""
    Z = differential_exp_basepoint(M, p, X, Y)
    differential_exp_basepoint!(M, Z, p, X, Y)

Compute ``D_p\exp_p X[Y]`` (in place of `Z`).

# See also
[`differential_exp_argument`](@ref), [`jacobi_field`](@ref)
"""
function differential_exp_basepoint(M::AbstractManifold, p, X, Y)
    return jacobi_field(M, p, exp(M, p, X), 1.0, Y, βdifferential_exp_basepoint)
end
function differential_exp_basepoint!(M::AbstractManifold, Z, p, X, Y)
    return jacobi_field!(M, Z, p, exp(M, p, X), 1.0, Y, βdifferential_exp_basepoint)
end

@doc raw"""
    Z = differential_exp_argument(M, p, X, Y)
    differential_exp_argument!(M, Z, p, X, Y)

computes ``D_X\exp_pX[Y]`` (in place of `Z`).
Note that ``X ∈  T_X(T_p\mathcal M) = T_p\mathcal M`` is still a tangent vector.

# See also
 [`differential_exp_basepoint`](@ref), [`jacobi_field`](@ref)
"""
function differential_exp_argument(M::AbstractManifold, p, X, Y)
    return jacobi_field(M, p, exp(M, p, X), 1.0, Y, βdifferential_exp_argument)
end
function differential_exp_argument!(M::AbstractManifold, Z, p, X, Y)
    return jacobi_field!(M, Z, p, exp(M, p, X), 1.0, Y, βdifferential_exp_argument)
end

@doc raw"""
    Y = differential_log_basepoint(M, p, q, X)
    differential_log_basepoint!(M, Y, p, q, X)

computes ``D_p\log_pq[X]`` (in place of `Y`).

# See also
 [`differential_log_argument`](@ref), [`jacobi_field`](@ref)
"""
function differential_log_basepoint(M::AbstractManifold, p, q, X)
    return jacobi_field(M, p, q, 0.0, X, βdifferential_log_basepoint)
end
function differential_log_basepoint!(M::AbstractManifold, Y, p, q, X)
    return jacobi_field!(M, Y, p, q, 0.0, X, βdifferential_log_basepoint)
end

@doc raw"""
    Y = differential_log_argument(M, p, q, X)
    differential_log_argument!(M, Y, p, q, X)

computes ``D_q\log_pq[X]`` (in place of `Y`).

# See also
 [`differential_log_basepoint`](@ref), [`jacobi_field`](@ref)
"""
function differential_log_argument(M::AbstractManifold, p, q, X)
    # order of p and q has to be reversed in this call, cf. Persch, 2018 Lemma 2.3
    return jacobi_field(M, q, p, 1.0, X, βdifferential_log_argument)
end
function differential_log_argument!(M::AbstractManifold, Y, p, q, X)
    # order of p and q has to be reversed in this call, cf. Persch, 2018 Lemma 2.3
    return jacobi_field!(M, Y, q, p, 1.0, X, βdifferential_log_argument)
end

@doc raw"""
    differential_exp_argument_lie_approx(M::AbstractManifold, p, X, Y; n)

Approximate differential of exponential map based on Lie group exponential. The formula
reads (see Theorem 1.7 of [^Helgason1978])
```math
D_X \exp_{p}(X)[Y] = (\mathrm{d}L_{\exp_e(X)})_e\left(\sum_{k=0}^{n}\frac{(-1)^k}{(k+1)!}(\operatorname{ad}_X)^k(Y)\right)
```
where ``(\operatorname{ad}_X)^k(Y)`` is defined recursively as ``(\operatorname{ad}_X)^0(Y) = Y``,
``\operatorname{ad}_X^{k+1}(Y) = [X, \operatorname{ad}_X^k(Y)]``.

[^Helgason1978]:
    > S. Helgason, Differential Geometry, Lie Groups, and Symmetric Spaces, First Edition.
    > Academic Press, 1978.
"""
function differential_exp_argument_lie_approx(M::AbstractManifold, p, X, Y; n = 20)
    Z = allocate(X)
    return differential_exp_argument_lie_approx!(M, Z, p, X, Y; n)
end

function differential_exp_argument_lie_approx! end

@doc raw"""
    differential_inverse_retract_argument_fd_approx(
        M::AbstractManifold,
        p,
        q,
        X;
        retr::AbstractRetractionMethod = default_retraction_method(M),
        invretr::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
        h::Real=sqrt(eps(eltype(X))),
    )

Approximate the differential of the inverse retraction `invretr` using a finite difference
formula (see Eq. (16) in [^Zimmermann2019]):
```math
\frac{\operatorname{retr}^{-1}_q(\operatorname{retr}_p(hX)) - \operatorname{retr}^{-1}_q(\operatorname{retr}_p(-hX))}{2h}
```
where ``h`` is the finite difference step `h`, ``\operatorname{retr}^{-1}`` is the inverse
retraction `invretr` and ``\operatorname{retr}`` is the retraction `retr`.

[^Zimmermann2019]:
    > R. Zimmermann, “Hermite interpolation and data processing errors on Riemannian matrix
    > manifolds,” arXiv:1908.05875 [cs, math], Sep. 2019,
    > Available: http://arxiv.org/abs/1908.05875
"""
function differential_inverse_retract_argument_fd_approx(
    M::AbstractManifold,
    p,
    q,
    X;
    retr::AbstractRetractionMethod = default_retraction_method(M),
    invretr::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
    h::Real = sqrt(eps(eltype(X))),
)
    Y = allocate_result(M, differential_inverse_retract_argument_fd_approx, X, p, q)
    differential_inverse_retract_argument_fd_approx!(M, Y, p, q, X; retr, invretr, h)
    return Y
end

function differential_inverse_retract_argument_fd_approx!(
    M::AbstractManifold,
    Y,
    p,
    q,
    X;
    retr::AbstractRetractionMethod = default_retraction_method(M),
    invretr::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),
    h::Real = sqrt(eps(eltype(X))),
)
    p_tmp = retract(M, q, h * X, retr)
    inverse_retract!(M, Y, p, p_tmp, invretr)
    retract!(M, p_tmp, q, -h * X, retr)
    X_tmp = inverse_retract(M, p, p_tmp, invretr)
    Y .-= X_tmp
    Y ./= 2 * h
    return Y
end

@doc raw"""
    differential_project_basepoint(M, p, X, Y)
    differential_project_basepoint!(M, Z, p, X, Y)

Compute the differential of the projection ``\operatorname{proj}_{T_p\mathcal M}(X)``

```math
D_p \operatorname{proj}_{T_p\mathcal M}(X)[Y],
```

with respect to the `basepoint` ``p``, i.e. ``X`` is a tangent vector in the embedding, ``Y \in T_p\mathcal M``, as well as the result.
"""
differential_project_basepoint(M::AbstractManifold, p, X, Y)

@doc raw"""
    differential_project_basepoint(M, p, X, Y)
    differential_project_basepoint!(M, Z, p, X, Y)

Compute the differential of the projection ``\operatorname{proj}_{T_p\mathcal M}(X)``
```math
D_p \operatorname{proj}_{T_p\mathcal M}(X)[Y]
```
with respect to the `basepoint` ``p`` in place of `Z`,
i.e. ``X`` is a tangent vector in the embedding, and ``Y,Z \in T_p\mathcal M``.
"""
differential_project_basepoint!(M::AbstractManifold, Z, p, X, Y)
