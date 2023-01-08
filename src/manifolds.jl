function adjoint_Jacobi_field(
    M::ProductManifold,
    p::ArrayPartition,
    q::ArrayPartition,
    t,
    X::ArrayPartition,
    β::Tβ,
) where {Tβ}
    return ArrayPartition(
        map(
            adjoint_Jacobi_field,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, q),
            ntuple(_ -> t, length(M.manifolds)),
            submanifold_components(M, X),
            ntuple(_ -> β, length(M.manifolds)),
        )...,
    )
end
function adjoint_Jacobi_field!(M::ProductManifold, Y, p, q, t, X, β::Tβ) where {Tβ}
    map(
        adjoint_Jacobi_field!,
        M.manifolds,
        submanifold_components(M, Y),
        submanifold_components(M, p),
        submanifold_components(M, q),
        ntuple(_ -> t, length(M.manifolds)),
        submanifold_components(M, X),
        ntuple(_ -> β, length(M.manifolds)),
    )
    return Y
end
function adjoint_Jacobi_field(::Circle{ℝ}, p, q, t, X, β::Tβ) where {Tβ}
    return X
end
function adjoint_Jacobi_field(::Euclidean{Tuple{}}, p, q, t, X, β::Tβ) where {Tβ}
    return X
end

function diagonalizing_projectors(M::AbstractSphere{ℝ}, p, X)
    X_norm = norm(M, p, X)
    X_normed = X / X_norm
    return (
        (zero(number_eltype(p)), ProjectorOntoVector(M, p, X_normed)),
        (one(number_eltype(p)), CoprojectorOntoVector(M, p, X_normed)),
    )
end

function diagonalizing_projectors(M::Hyperbolic, p, X)
    X_norm = norm(M, p, X)
    X_normed = X / X_norm
    return (
        (zero(number_eltype(p)), ProjectorOntoVector(M, p, X_normed)),
        (-one(number_eltype(p)), CoprojectorOntoVector(M, p, X_normed)),
    )
end

function diagonalizing_projectors(::Euclidean, p, X)
    return ((zero(number_eltype(p)), IdentityProjector()),)
end

function diagonalizing_projectors(M::Circle{ℝ}, p, X)
    sbv = sign(X[])
    proj = ProjectorOntoVector(
        M,
        p,
        Manifolds.StaticArrays.@SVector [sbv == 0 ? one(sbv) : sbv]
    )
    return ((zero(number_eltype(p)), proj),)
end

function jacobi_field(
    M::ProductManifold,
    p::ArrayPartition,
    q::ArrayPartition,
    t,
    X::ArrayPartition,
    β::Tβ,
) where {Tβ}
    return ArrayPartition(
        map(
            jacobi_field,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, q),
            ntuple(_ -> t, length(M.manifolds)),
            submanifold_components(M, X),
            ntuple(_ -> β, length(M.manifolds)),
        )...,
    )
end
function jacobi_field!(M::ProductManifold, Y, p, q, t, X, β::Tβ) where {Tβ}
    map(
        jacobi_field!,
        M.manifolds,
        submanifold_components(M, Y),
        submanifold_components(M, p),
        submanifold_components(M, q),
        ntuple(_ -> t, length(M.manifolds)),
        submanifold_components(M, X),
        ntuple(_ -> β, length(M.manifolds)),
    )
    return Y
end
function jacobi_field(::Circle{ℝ}, p, q, t, X, β::Tβ) where {Tβ}
    return X
end
function jacobi_field(::Euclidean{Tuple{}}, p, q, t, X, β::Tβ) where {Tβ}
    return X
end
