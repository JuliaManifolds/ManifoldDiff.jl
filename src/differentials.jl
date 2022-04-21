
@doc raw"""
    differential_geodesic_startpoint(M, p, q, t, X)
    differential_geodesic_startpoint!(M, Y, p, q, t, X)

Compute ``D_p g(t;p,q)[η]`` (in place of `Y`).

# See also
[`differential_geodesic_endpoint`](@ref), [`jacobi_field`](@ref)
"""
function differential_geodesic_startpoint(M::AbstractManifold, p, q, t, X)
    return jacobi_field(M, p, q, t, X, βdifferential_geodesic_startpoint)
end
function differential_geodesic_startpoint!(M::AbstractManifold, Y, p, q, t, X)
    return jacobi_field!(M, Y, p, q, t, X, βdifferential_geodesic_startpoint)
end


@doc raw"""
    differential_geodesic_endpoint(M, p, q, t, X)
    differential_geodesic_endpoint!(M, Y, p, q, t, X)

Compute ``D_qg(t;p,q)[X]`` (in place of `Y`).

# See also
 [`differential_geodesic_startpoint`](@ref), [`jacobi_field`](@ref)
"""
function differential_geodesic_endpoint(M::AbstractManifold, p, q, t, X)
    return jacobi_field(M, q, p, 1 - t, X, βdifferential_geodesic_startpoint)
end
function differential_geodesic_endpoint!(M::AbstractManifold, Y, p, q, t, X)
    return jacobi_field!(M, Y, q, p, 1 - t, X, βdifferential_geodesic_startpoint)
end

@doc raw"""
    differential_exp_basepoint(M, p, X, Y)
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
    differential_exp_argument(M, p, X, Y)
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
    differential_log_basepoint(M, p, q, X)
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
    differential_log_argument(M, p, q, X)
    differential_log_argument(M, Y, p, q, X)

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
