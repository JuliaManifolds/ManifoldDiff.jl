
@doc raw"""
    abstract type AbstractProjector end

An abstract type for projectors on a tangent space $T_pM$ for fixed values of `p` and `M`.
Calling a projector on a tangent vector returns a new tangent vector:

    (Π::AbstractProjector)(X) -> Y

Projectors assume that `X` is a valid vector from $T_pM$.
"""
abstract type AbstractProjector end

@doc raw"""
    diagonalizing_projectors(M::AbstractManifold, p, X)

Compute eigenvalues of the Jacobi operator $Y → R(Y,X)X$, where $R$ is the curvature
endomorphism, together with projectors onto eigenspaces of the operator.
Projectors are objects of subtypes of [`AbstractProjector`](@ref).

By default constructs projectors using the [`DiagonalizingOrthonormalBasis`](@ref).
"""
function diagonalizing_projectors(M::AbstractManifold, p, X)
    B = get_basis(M, p, DiagonalizingOrthonormalBasis(X))
    N = length(B.data.eigenvalues)
    return map(
        i -> (B.data.eigenvalues[i], ProjectorOntoVector(M, p, B.data.vectors[i])),
        1:N,
    )
end

"""
    ProjectorOntoVector{TM<:AbstractManifold,TP,TX}
    
A structure that represents projector onto the subspace of the tangent space at `p`
from manifold `M` spanned by tangent vector `X` of unit norm.

# Constructor

ProjectorOntoVector(M::AbstractManifold, p, X)
"""
struct ProjectorOntoVector{TM<:AbstractManifold,TP,TX} <: AbstractProjector
    M::TM
    p::TP
    X::TX
end

function (Π::ProjectorOntoVector)(X)
    return inner(Π.M, Π.p, X, Π.X) * Π.X
end

"""
    CoprojectorOntoVector{TM<:AbstractManifold,TP,TX}
    

A structure that represents projector onto the subspace of the tangent space at `p`
from manifold `M` othogonal to vector `X` of unit norm.

# Constructor

CoprojectorOntoVector(M::AbstractManifold, p, X)
"""
struct CoprojectorOntoVector{TM<:AbstractManifold,TP,TX} <: AbstractProjector
    M::TM
    p::TP
    X::TX
end

function (Π::CoprojectorOntoVector)(X)
    return X - inner(Π.M, Π.p, X, Π.X) * Π.X
end

struct IdentityProjector <: AbstractProjector end

(::IdentityProjector)(X) = X
