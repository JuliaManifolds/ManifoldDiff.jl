
function allocate_jacobian(
    M::AbstractManifold,
    f,
    p;
    Barg::AbstractBasis = DefaultOrthonormalBasis(),
    Bval::AbstractBasis = DefaultOrthonormalBasis(),
)
    n = number_of_coordinates(M, Barg)
    m = number_of_coordinates(M, Bval)
    J = zeros(float(number_eltype(p)), n, m)
    return J
end

function jacobian_exp_argument!(
    M::AbstractManifold,
    J::AbstractMatrix,
    p,
    X,
    Barg::AbstractBasis = DefaultOrthonormalBasis(),
    Bval::AbstractBasis = DefaultOrthonormalBasis(),
)
    Bb = get_basis(M, p, Barg)
    Vs = get_vectors(M, p, Bb)
    Z = similar(X)
    p_exp = exp(M, p, X)
    for i in 1:length(Vs)
        differential_exp_argument!(M, Z, p, X, Vs[i])
        get_coordinates!(M, view(J, :, i), p_exp, Z, Bval)
    end
    return J
end

@doc raw"""
    jacobian_exp_argument(
        M::AbstractManifold,
        p,
        X,
        Barg::AbstractBasis=DefaultOrthonormalBasis(),
        Bval::AbstractBasis=DefaultOrthonormalBasis(),
    )

Compute Jacobian of the exponential map with respect to its argument (tangent vector).
Differential of the exponential map is here considered as a function from ``T_p \mathcal{M}``
to ``T_{\exp_p X} \mathcal{M}``. Jacobian coefficients are represented in basis `Barg`
in the domain and in `Bval` in the codomain.
"""
function jacobian_exp_argument(
    M::AbstractManifold,
    p,
    X;
    Barg::AbstractBasis = DefaultOrthonormalBasis(),
    Bval::AbstractBasis = DefaultOrthonormalBasis(),
)
    J = allocate_jacobian(M, jacobian_exp_argument, X; Barg = Barg, Bval = Bval)
    jacobian_exp_argument!(M, J, p, X, Barg, Bval)
    return J
end

function jacobian_exp_basepoint!(
    M::AbstractManifold,
    J::AbstractMatrix,
    p,
    X,
    Barg::AbstractBasis = DefaultOrthonormalBasis(),
    Bval::AbstractBasis = DefaultOrthonormalBasis(),
)
    Bb = get_basis(M, p, Barg)
    Vs = get_vectors(M, p, Bb)
    Z = similar(X)
    p_exp = exp(M, p, X)
    for i in 1:length(Vs)
        differential_exp_basepoint!(M, Z, p, X, Vs[i])
        get_coordinates!(M, view(J, :, i), p_exp, Z, Bval)
    end
    return J
end

@doc raw"""
    jacobian_exp_basepoint(
        M::AbstractManifold,
        p,
        X,
        Barg::AbstractBasis=DefaultOrthonormalBasis(),
        Bval::AbstractBasis=DefaultOrthonormalBasis(),
    )

Compute Jacobian of the exponential map with respect to the basepoint.
Differential of the exponential map is here considered as a function from ``T_p \mathcal{M}``
to ``T_{\exp_p X} \mathcal{M}``. Jacobian coefficients are represented in basis `Barg`
in the domain and in `Bval` in the codomain.
"""
function jacobian_exp_basepoint(
    M::AbstractManifold,
    p,
    X,
    Barg::AbstractBasis = DefaultOrthonormalBasis(),
    Bval::AbstractBasis = DefaultOrthonormalBasis(),
)
    J = allocate_jacobian(M, jacobian_exp_basepoint, X; Barg = Barg, Bval = Bval)
    jacobian_exp_basepoint!(M, J, p, X, Barg, Bval)
    return J
end

function jacobian_log_argument!(
    M::AbstractManifold,
    J::AbstractMatrix,
    p,
    q,
    Barg::AbstractBasis = DefaultOrthonormalBasis(),
    Bval::AbstractBasis = DefaultOrthonormalBasis(),
)
    Bb = get_basis(M, q, Barg)
    Vs = get_vectors(M, q, Bb)
    Z = zero_vector(M, p)
    for i in 1:length(Vs)
        differential_log_argument!(M, Z, p, q, Vs[i])
        get_coordinates!(M, view(J, :, i), p, Z, Bval)
    end
    return J
end

@doc raw"""
    jacobian_log_argument(
        M::AbstractManifold,
        p,
        q,
        Barg::AbstractBasis=DefaultOrthonormalBasis(),
        Bval::AbstractBasis=DefaultOrthonormalBasis(),
    )

Compute Jacobian of the logarithmic map with respect to its argument (point `q`).
Differential of the logarithmic map is here considered as a function from ``T_q \mathcal{M}``
to ``T_p \mathcal{M}``. Jacobian coefficients are represented in basis `Barg`
in the domain and in `Bval` in the codomain.
"""
function jacobian_log_argument(
    M::AbstractManifold,
    p,
    q;
    Barg::AbstractBasis = DefaultOrthonormalBasis(),
    Bval::AbstractBasis = DefaultOrthonormalBasis(),
)
    J = allocate_jacobian(M, jacobian_log_argument, q; Barg = Barg, Bval = Bval)
    jacobian_log_argument!(M, J, p, q, Barg, Bval)
    return J
end

function jacobian_log_basepoint!(
    M::AbstractManifold,
    J::AbstractMatrix,
    p,
    q,
    Barg::AbstractBasis = DefaultOrthonormalBasis(),
    Bval::AbstractBasis = DefaultOrthonormalBasis(),
)
    Bb = get_basis(M, p, Barg)
    Vs = get_vectors(M, p, Bb)
    Z = zero_vector(M, p)
    for i in 1:length(Vs)
        differential_log_basepoint!(M, Z, p, q, Vs[i])
        get_coordinates!(M, view(J, :, i), p, Z, Bval)
    end
    return J
end

@doc raw"""
    jacobian_log_basepoint(
        M::AbstractManifold,
        p,
        q,
        Barg::AbstractBasis=DefaultOrthonormalBasis(),
        Bval::AbstractBasis=DefaultOrthonormalBasis(),
    )

Compute Jacobian of the logarithmic map with respect to the basepoint.
Differential of the logarithmic map is here considered as a function from ``T_q \mathcal{M}``
to ``T_p \mathcal{M}``. Jacobian coefficients are represented in basis `Barg`
in the domain and in `Bval` in the codomain.
"""
function jacobian_log_basepoint(
    M::AbstractManifold,
    p,
    q,
    Barg::AbstractBasis = DefaultOrthonormalBasis(),
    Bval::AbstractBasis = DefaultOrthonormalBasis(),
)
    J = allocate_jacobian(M, jacobian_log_basepoint, q; Barg = Barg, Bval = Bval)
    jacobian_log_basepoint!(M, J, p, q, Barg, Bval)
    return J
end
