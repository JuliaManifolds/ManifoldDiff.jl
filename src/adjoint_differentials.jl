
@doc raw"""
    adjoint_differential_geodesic_startpoint(M, p, q, t, X)
    adjoint_differential_geodesic_startpoint!(M, Y, p, q, t, X)

Compute the adjoint of ``D_p γ(t; p, q)[X]`` (in place of `Y`).

# See also

[`differential_geodesic_startpoint`](@ref), [`adjoint_Jacobi_field`](@ref)
"""
function adjoint_differential_geodesic_startpoint(M::AbstractManifold, p, q, t, X)
    return adjoint_Jacobi_field(M, p, q, t, X, βdifferential_geodesic_startpoint)
end
function adjoint_differential_geodesic_startpoint!(M::AbstractManifold, Y, p, q, t, X)
    return adjoint_Jacobi_field!(M, Y, p, q, t, X, βdifferential_geodesic_startpoint)
end

@doc raw"""
    adjoint_differential_geodesic_endpoint(M, p, q, t, X)
    adjoint_differential_geodesic_endpoint!(M, Y, p, q, t, X)

Compute the adjoint of ``D_q γ(t; p, q)[X]`` (in place of `Y`).

# See also

[`differential_geodesic_endpoint`](@ref), [`adjoint_Jacobi_field`](@ref)
"""
function adjoint_differential_geodesic_endpoint(M::AbstractManifold, p, q, t, X)
    return adjoint_Jacobi_field(M, q, p, 1 - t, X, βdifferential_geodesic_startpoint)
end
function adjoint_differential_geodesic_endpoint!(M::AbstractManifold, Y, p, q, t, X)
    return adjoint_Jacobi_field!(M, Y, q, p, 1 - t, X, βdifferential_geodesic_startpoint)
end

@doc raw"""
    adjoint_differential_exp_basepoint(M, p, X, Y)
    adjoint_differential_exp_basepoint!(M, Z, p, X, Y)

Computes the adjoint of ``D_p \exp_p X[Y]`` (in place of `Z`).

# See also

[`differential_exp_basepoint`](@ref), [`adjoint_Jacobi_field`](@ref)
"""
function adjoint_differential_exp_basepoint(M::AbstractManifold, p, X, Y)
    return adjoint_Jacobi_field(M, p, exp(M, p, X), 1.0, Y, βdifferential_exp_basepoint)
end
function adjoint_differential_exp_basepoint!(M::AbstractManifold, Z, p, X, Y)
    return adjoint_Jacobi_field!(M, Z, p, exp(M, p, X), 1.0, Y, βdifferential_exp_basepoint)
end

@doc raw"""
    adjoint_differential_exp_argument(M, p, X, Y)
    adjoint_differential_exp_argument!(M, Z, p, X, Y)

Compute the adjoint of ``D_X\exp_p X[Y]`` (in place of `Z`).
Note that ``X ∈  T_p(T_p\mathcal M) = T_p\mathcal M`` is still a tangent vector.

# See also

[`differential_exp_argument`](@ref), [`adjoint_Jacobi_field`](@ref)
"""
function adjoint_differential_exp_argument(M::AbstractManifold, p, X, Y)
    return adjoint_Jacobi_field(M, p, exp(M, p, X), 1.0, Y, βdifferential_exp_argument)
end
function adjoint_differential_exp_argument!(M::AbstractManifold, Z, p, X, Y)
    return adjoint_Jacobi_field!(M, Z, p, exp(M, p, X), 1.0, Y, βdifferential_exp_argument)
end

@doc raw"""
    adjoint_differential_log_basepoint(M, p, q, X)
    adjoint_differential_log_basepoint!(M, Y, p, q, X)

computes the adjoint of ``D_p log_p q[X]`` (in place of `Y`).

# See also
[`differential_log_basepoint`](@ref), [`adjoint_Jacobi_field`](@ref)
"""
function adjoint_differential_log_basepoint(M::AbstractManifold, p, q, X)
    return adjoint_Jacobi_field(M, p, q, 0.0, X, βdifferential_log_basepoint)
end
function adjoint_differential_log_basepoint!(M::AbstractManifold, Y, p, q, X)
    return adjoint_Jacobi_field!(M, Y, p, q, 0.0, X, βdifferential_log_basepoint)
end

@doc raw"""
    adjoint_differential_log_argument(M, p, q, X)
    adjoint_differential_log_argument!(M, Y, p, q, X)

Compute the adjoint of ``D_q log_p q[X]`` (in place of `Y`).

# See also
[`differential_log_argument`](@ref), [`adjoint_Jacobi_field`](@ref)
"""
function adjoint_differential_log_argument(M::AbstractManifold, p, q, X)
    # order of p and q has to be reversed in this call, cf. Persch, 2018 Lemma 2.3
    return adjoint_Jacobi_field(M, q, p, 1.0, X, βdifferential_log_argument)
end
function adjoint_differential_log_argument!(M::AbstractManifold, Y, p, q, X)
    return adjoint_Jacobi_field!(M, Y, q, p, 1.0, X, βdifferential_log_argument)
end
