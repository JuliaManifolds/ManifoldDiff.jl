module ManifoldDiff

using Markdown: @doc_str
using LinearAlgebra

using StaticArrays

using Manifolds
using Manifolds:
    AbstractDiffBackend,
    NoneDiffBackend,
    default_differential_backend,
    set_default_differential_backend!
using ManifoldsBase
using ChainRulesCore
using ForwardDiff
using RecursiveArrayTools
using ManifoldsBase:
    TangentSpaceType,
    PowerManifoldNested,
    PowerManifoldNestedReplacing,
    get_iterator,
    _write,
    _read


using Requires

import Manifolds:
    _derivative, _derivative!, _gradient, _gradient!, _jacobian, _jacobian!, _hessian

"""
    hessian(M::AbstractManifold, f, p, backend::Manifolds.TangentDiffBackend)

Compute the Hessian of function `f` at point `p` using the given `backend`. The formula
for normal coordinate systems from[^SommerFletcherPennec2020] is used.

[^SommerFletcherPennec2020]:
    > S. Sommer, T. Fletcher, and X. Pennec, “1 - Introduction to differential and Riemannian
    > geometry,” in Riemannian Geometric Statistics in Medical Image Analysis, X. Pennec,
    > S. Sommer, and T. Fletcher, Eds. Academic Press, 2020, pp. 3–37.
    > doi: 10.1016/B978-0-12-814725-2.00008-X.
"""
function hessian(M::AbstractManifold, f, p, backend::Manifolds.TangentDiffBackend)
    X = get_coordinates(M, p, zero_vector(M, p), backend.basis_arg)
    onb_coords = _hessian(X, backend.diff_backend) do Y
        return f(exp(M, p, get_vector(M, p, Y, backend.basis_arg)))
    end
    return onb_coords
end


include("diagonalizing_projectors.jl")

include("differentials.jl")
include("adjoint_differentials.jl")
include("Jacobi_fields.jl")

function __init__()
    @require FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41" begin
        using .FiniteDiff
        include("finite_diff.jl")
    end

    @require ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210" begin
        using .ForwardDiff
        include("forward_diff.jl")
    end

    @require ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267" begin
        using .ReverseDiff: ReverseDiff
        include("reverse_diff.jl")
    end

    @require Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f" begin
        using .Zygote: Zygote
        include("zygote.jl")
    end

    return nothing
end


end # module
