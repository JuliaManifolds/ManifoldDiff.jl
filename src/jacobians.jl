
"""
    allocate_jacobian(
        M_domain::AbstractManifold,
        M_codomain::AbstractManifold,
        f,
        p;
        basis_domain::AbstractBasis = DefaultOrthonormalBasis(),
        basis_codomain::AbstractBasis = DefaultOrthonormalBasis(),
    )

Allocate Jacobian of function `f` with given domain and codomain at point `p`.
`basis_domain` and `basis_codomain` denote bases of tangent spaces at, respectively, `p` and `f(p)`.
"""
function allocate_jacobian(
    M_domain::AbstractManifold,
    M_codomain::AbstractManifold,
    f,
    p;
    basis_domain::AbstractBasis = DefaultOrthonormalBasis(),
    basis_codomain::AbstractBasis = DefaultOrthonormalBasis(),
)
    n = number_of_coordinates(M_domain, basis_domain)
    m = number_of_coordinates(M_codomain, basis_codomain)
    J = zeros(float(number_eltype(p)), n, m)
    return J
end

function jacobian_exp_argument!(
    M::AbstractManifold,
    J::AbstractMatrix,
    p,
    X,
    basis_domain::AbstractBasis = DefaultOrthonormalBasis(),
    basis_codomain::AbstractBasis = DefaultOrthonormalBasis(),
)
    Bb = get_basis(M, p, basis_domain)
    Vs = get_vectors(M, p, Bb)
    Z = similar(X)
    p_exp = exp(M, p, X)
    for i in 1:length(Vs)
        differential_exp_argument!(M, Z, p, X, Vs[i])
        get_coordinates!(M, view(J, :, i), p_exp, Z, basis_codomain)
    end
    return J
end

@doc raw"""
    jacobian_exp_argument(
        M::AbstractManifold,
        p,
        X,
        basis_domain::AbstractBasis=DefaultOrthonormalBasis(),
        basis_codomain::AbstractBasis=DefaultOrthonormalBasis(),
    )

Compute Jacobian of the exponential map with respect to its argument (tangent vector).
Differential of the exponential map is here considered as a function from ``T_p \mathcal{M}``
to ``T_{\exp_p X} \mathcal{M}``. Jacobian coefficients are represented in basis `basis_domain`
in the domain and in `basis_codomain` in the codomain.
"""
function jacobian_exp_argument(
    M::AbstractManifold,
    p,
    X,
    basis_domain::AbstractBasis = DefaultOrthonormalBasis(),
    basis_codomain::AbstractBasis = DefaultOrthonormalBasis(),
)
    J = allocate_jacobian(
        M,
        M,
        jacobian_exp_argument,
        X;
        basis_domain = basis_domain,
        basis_codomain = basis_codomain,
    )
    jacobian_exp_argument!(M, J, p, X, basis_domain, basis_codomain)
    return J
end

function jacobian_exp_basepoint!(
    M::AbstractManifold,
    J::AbstractMatrix,
    p,
    X,
    basis_domain::AbstractBasis = DefaultOrthonormalBasis(),
    basis_codomain::AbstractBasis = DefaultOrthonormalBasis(),
)
    Bb = get_basis(M, p, basis_domain)
    Vs = get_vectors(M, p, Bb)
    Z = similar(X)
    p_exp = exp(M, p, X)
    for i in 1:length(Vs)
        differential_exp_basepoint!(M, Z, p, X, Vs[i])
        get_coordinates!(M, view(J, :, i), p_exp, Z, basis_codomain)
    end
    return J
end

@doc raw"""
    jacobian_exp_basepoint(
        M::AbstractManifold,
        p,
        X,
        basis_domain::AbstractBasis=DefaultOrthonormalBasis(),
        basis_codomain::AbstractBasis=DefaultOrthonormalBasis(),
    )

Compute Jacobian of the exponential map with respect to the basepoint.
Differential of the exponential map is here considered as a function from ``T_p \mathcal{M}``
to ``T_{\exp_p X} \mathcal{M}``. Jacobian coefficients are represented in basis `basis_domain`
in the domain and in `basis_codomain` in the codomain.
"""
function jacobian_exp_basepoint(
    M::AbstractManifold,
    p,
    X,
    basis_domain::AbstractBasis = DefaultOrthonormalBasis(),
    basis_codomain::AbstractBasis = DefaultOrthonormalBasis(),
)
    J = allocate_jacobian(
        M,
        M,
        jacobian_exp_basepoint,
        X;
        basis_domain = basis_domain,
        basis_codomain = basis_codomain,
    )
    jacobian_exp_basepoint!(M, J, p, X, basis_domain, basis_codomain)
    return J
end

function jacobian_log_argument!(
    M::AbstractManifold,
    J::AbstractMatrix,
    p,
    q,
    basis_domain::AbstractBasis = DefaultOrthonormalBasis(),
    basis_codomain::AbstractBasis = DefaultOrthonormalBasis(),
)
    Bb = get_basis(M, q, basis_domain)
    Vs = get_vectors(M, q, Bb)
    Z = zero_vector(M, p)
    for i in 1:length(Vs)
        differential_log_argument!(M, Z, p, q, Vs[i])
        get_coordinates!(M, view(J, :, i), p, Z, basis_codomain)
    end
    return J
end

@doc raw"""
    jacobian_log_argument(
        M::AbstractManifold,
        p,
        q,
        basis_domain::AbstractBasis=DefaultOrthonormalBasis(),
        basis_codomain::AbstractBasis=DefaultOrthonormalBasis(),
    )

Compute Jacobian of the logarithmic map with respect to its argument (point `q`).
Differential of the logarithmic map is here considered as a function from ``T_q \mathcal{M}``
to ``T_p \mathcal{M}``. Jacobian coefficients are represented in basis `basis_domain`
in the domain and in `basis_codomain` in the codomain.
"""
function jacobian_log_argument(
    M::AbstractManifold,
    p,
    q,
    basis_domain::AbstractBasis = DefaultOrthonormalBasis(),
    basis_codomain::AbstractBasis = DefaultOrthonormalBasis(),
)
    J = allocate_jacobian(
        M,
        M,
        jacobian_log_argument,
        q;
        basis_domain = basis_domain,
        basis_codomain = basis_codomain,
    )
    jacobian_log_argument!(M, J, p, q, basis_domain, basis_codomain)
    return J
end

function jacobian_log_basepoint!(
    M::AbstractManifold,
    J::AbstractMatrix,
    p,
    q,
    basis_domain::AbstractBasis = DefaultOrthonormalBasis(),
    basis_codomain::AbstractBasis = DefaultOrthonormalBasis(),
)
    Bb = get_basis(M, p, basis_domain)
    Vs = get_vectors(M, p, Bb)
    Z = zero_vector(M, p)
    for i in 1:length(Vs)
        differential_log_basepoint!(M, Z, p, q, Vs[i])
        get_coordinates!(M, view(J, :, i), p, Z, basis_codomain)
    end
    return J
end

@doc raw"""
    jacobian_log_basepoint(
        M::AbstractManifold,
        p,
        q,
        basis_domain::AbstractBasis=DefaultOrthonormalBasis(),
        basis_codomain::AbstractBasis=DefaultOrthonormalBasis(),
    )

Compute Jacobian of the logarithmic map with respect to the basepoint.
Differential of the logarithmic map is here considered as a function from ``T_q \mathcal{M}``
to ``T_p \mathcal{M}``. Jacobian coefficients are represented in basis `basis_domain`
in the domain and in `basis_codomain` in the codomain.
"""
function jacobian_log_basepoint(
    M::AbstractManifold,
    p,
    q,
    basis_domain::AbstractBasis = DefaultOrthonormalBasis(),
    basis_codomain::AbstractBasis = DefaultOrthonormalBasis(),
)
    J = allocate_jacobian(
        M,
        M,
        jacobian_log_basepoint,
        q;
        basis_domain = basis_domain,
        basis_codomain = basis_codomain,
    )
    jacobian_log_basepoint!(M, J, p, q, basis_domain, basis_codomain)
    return J
end
