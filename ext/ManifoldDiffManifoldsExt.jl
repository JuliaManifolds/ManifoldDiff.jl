module ManifoldDiffManifoldsExt

if isdefined(Base, :get_extension)
    using ManifoldDiff
    using ManifoldDiff: ProjectorOntoVector, CoprojectorOntoVector, IdentityProjector
    using Manifolds
else
    # imports need to be relative for Requires.jl-based workflows:
    # https://github.com/JuliaArrays/ArrayInterface.jl/pull/387
    using ..ManifoldDiff
    using ..ManifoldDiff: ProjectorOntoVector, CoprojectorOntoVector, IdentityProjector
    using ..Manifolds
end

function ManifoldDiff.adjoint_Jacobi_field(
    M::ProductManifold,
    p::ArrayPartition,
    q::ArrayPartition,
    t,
    X::ArrayPartition,
    β::Tβ,
) where {Tβ}
    return ArrayPartition(
        map(
            ManifoldDiff.adjoint_Jacobi_field,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, q),
            ntuple(_ -> t, length(M.manifolds)),
            submanifold_components(M, X),
            ntuple(_ -> β, length(M.manifolds)),
        )...,
    )
end
function ManifoldDiff.adjoint_Jacobi_field!(
    M::ProductManifold,
    Y,
    p,
    q,
    t,
    X,
    β::Tβ,
) where {Tβ}
    map(
        ManifoldDiff.adjoint_Jacobi_field!,
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
function ManifoldDiff.adjoint_Jacobi_field(::Circle{ℝ}, p, q, t, X, β::Tβ) where {Tβ}
    return X
end
function ManifoldDiff.adjoint_Jacobi_field(
    ::Euclidean{Tuple{}},
    p,
    q,
    t,
    X,
    β::Tβ,
) where {Tβ}
    return X
end

function ManifoldDiff.diagonalizing_projectors(M::AbstractSphere{ℝ}, p, X)
    X_norm = norm(M, p, X)
    X_normed = X / X_norm
    return (
        (zero(number_eltype(p)), ProjectorOntoVector(M, p, X_normed)),
        (one(number_eltype(p)), CoprojectorOntoVector(M, p, X_normed)),
    )
end

function ManifoldDiff.diagonalizing_projectors(M::Hyperbolic, p, X)
    X_norm = norm(M, p, X)
    X_normed = X / X_norm
    return (
        (zero(number_eltype(p)), ProjectorOntoVector(M, p, X_normed)),
        (-one(number_eltype(p)), CoprojectorOntoVector(M, p, X_normed)),
    )
end

function ManifoldDiff.diagonalizing_projectors(::Euclidean, p, X)
    return ((zero(number_eltype(p)), IdentityProjector()),)
end

function ManifoldDiff.diagonalizing_projectors(M::Circle{ℝ}, p, X)
    sbv = sign(X[])
    proj = ProjectorOntoVector(
        M,
        p,
        Manifolds.StaticArrays.@SVector [sbv == 0 ? one(sbv) : sbv]
    )
    return ((zero(number_eltype(p)), proj),)
end

function ManifoldDiff.jacobi_field(
    M::ProductManifold,
    p::ArrayPartition,
    q::ArrayPartition,
    t,
    X::ArrayPartition,
    β::Tβ,
) where {Tβ}
    return ArrayPartition(
        map(
            ManifoldDiff.jacobi_field,
            M.manifolds,
            submanifold_components(M, p),
            submanifold_components(M, q),
            ntuple(_ -> t, length(M.manifolds)),
            submanifold_components(M, X),
            ntuple(_ -> β, length(M.manifolds)),
        )...,
    )
end
function ManifoldDiff.jacobi_field!(M::ProductManifold, Y, p, q, t, X, β::Tβ) where {Tβ}
    map(
        ManifoldDiff.jacobi_field!,
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
function ManifoldDiff.jacobi_field(::Circle{ℝ}, p, q, t, X, β::Tβ) where {Tβ}
    return X
end
function ManifoldDiff.jacobi_field(::Euclidean{Tuple{}}, p, q, t, X, β::Tβ) where {Tβ}
    return X
end

### differentials

function ManifoldDiff.differential_exp_argument_lie_approx!(
    M::AbstractManifold,
    Z,
    p,
    X,
    Y;
    n = 20,
)
    tmp = copy(M, p, Y)
    a = -1.0
    zero_vector!(M, Z, p)
    for k in 0:n
        a *= -1 // (k + 1)
        Z .+= a .* tmp
        if k < n
            copyto!(tmp, lie_bracket(M, X, tmp))
        end
    end
    q = exp(M, p, X)
    translate_diff!(M, Z, q, Identity(M), Z)
    return Z
end


end
