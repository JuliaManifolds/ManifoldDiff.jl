module ManifoldDiff

using Markdown: @doc_str
using LinearAlgebra

using ManifoldsBase
# using ChainRulesCore
# using RecursiveArrayTools

using ManifoldsBase:
    AbstractBasis,
    AbstractMetric,
    EuclideanMetric,
    TangentSpaceType,
    PowerManifoldNested,
    PowerManifoldNestedReplacing,
    get_iterator,
    _write,
    _read


using Requires


"""
    AbstractDiffBackend

An abstract type for diff backends. See [`FiniteDifferencesBackend`](@ref) for
an example.
"""
abstract type AbstractDiffBackend end

struct NoneDiffBackend <: AbstractDiffBackend end

"""
    _derivative(f, t[, backend::AbstractDiffBackend])

Compute the derivative of a callable `f` at time `t` computed using the given `backend`,
an object of type [`AbstractDiffBackend`](@ref ManifoldDiff.AbstractDiffBackend). If the backend is not explicitly
specified, it is obtained using the function [`default_differential_backend`](@ref).

This function calculates plain Euclidean derivatives, for Riemannian differentiation see
for example [`differential`](@ref differential(::AbstractManifold, ::Any, ::Real, ::AbstractRiemannianDiffBackend)).

!!! note

    Not specifying the backend explicitly will usually result in a type instability
    and decreased performance.
"""
function _derivative end

_derivative(f, t) = _derivative(f, t, default_differential_backend())

function _derivative!(
    f,
    X,
    t,
    backend::AbstractDiffBackend = default_differential_backend(),
)
    return copyto!(X, _derivative(f, t, backend))
end

"""
    _gradient(f, p[, backend::AbstractDiffBackend])

Compute the gradient of a callable `f` at point `p` computed using the given `backend`,
an object of type [`AbstractDiffBackend`](@ref ManifoldDiff.AbstractDiffBackend). If the backend is not explicitly
specified, it is obtained using the function [`default_differential_backend`](@ref).

This function calculates plain Euclidean gradients, for Riemannian gradient calculation see
for example [`gradient`](@ref gradient(::AbstractManifold, ::Any, ::Any, ::AbstractRiemannianDiffBackend)).

!!! note

    Not specifying the backend explicitly will usually result in a type instability
    and decreased performance.
"""
function _gradient end

_gradient(f, p) = _gradient(f, p, default_differential_backend())

function _gradient!(f, X, p, backend::AbstractDiffBackend = default_differential_backend())
    return copyto!(X, _gradient(f, p, backend))
end

"""
    _hessian(f, p[, backend::AbstractDiffBackend])

Compute the Hessian of a callable `f` at point `p` computed using the given `backend`,
an object of type [`AbstractDiffBackend`](@ref). If the backend is not explicitly
specified, it is obtained using the function [`default_differential_backend`](@ref).

This function calculates plain Euclidean Hessian.

!!! note

    Not specifying the backend explicitly will usually result in a type instability
    and decreased performance.
"""
function _hessian end

_hessian(f, p) = _hessian(f, p, default_differential_backend())

"""
    _jacobian(f, p[, backend::AbstractDiffBackend])

Compute the Jacobian of a callable `f` at point `p` computed using the given `backend`,
an object of type [`AbstractDiffBackend`](@ref). If the backend is not explicitly
specified, it is obtained using the function [`default_differential_backend`](@ref).

This function calculates plain Euclidean Jacobians, for Riemannian Jacobian calculation see
for example [`gradient`](@ref gradient(::AbstractManifold, ::Any, ::Any, ::AbstractRiemannianDiffBackend)).

!!! note

    Not specifying the backend explicitly will usually result in a type instability
    and decreased performance.
"""
function _jacobian end

_jacobian(f, p) = _jacobian(f, p, default_differential_backend())

function _jacobian!(f, X, p, backend::AbstractDiffBackend = default_differential_backend())
    return copyto!(X, _jacobian(f, p, backend))
end

"""
    CurrentDiffBackend(backend::AbstractDiffBackend)

A mutable struct for storing the current differentiation backend in a global
constant [`_current_default_differential_backend`](@ref).

# See also

[`AbstractDiffBackend`](@ref), [`default_differential_backend`](@ref), [`set_default_differential_backend!`](@ref)
"""
mutable struct CurrentDiffBackend
    backend::AbstractDiffBackend
end

"""
    _current_default_differential_backend

The instance of [`CurrentDiffBackend`](@ref) that stores the globally default
differentiation backend.
"""
const _current_default_differential_backend = CurrentDiffBackend(NoneDiffBackend())
"""
    default_differential_backend() -> AbstractDiffBackend

Get the default differentiation backend.
"""
default_differential_backend() = _current_default_differential_backend.backend

"""
    set_default_differential_backend!(backend::AbstractDiffBackend)

Set current backend for differentiation to `backend`.
"""
function set_default_differential_backend!(backend::AbstractDiffBackend)
    _current_default_differential_backend.backend = backend
    return backend
end

include("diagonalizing_projectors.jl")

include("adjoint_differentials.jl")
include("derivatives.jl")
include("differentials.jl")
include("gradients.jl")
include("Jacobi_fields.jl")

include("riemannian_diff.jl")
include("embedded_diff.jl")


function __init__()
    @require FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41" begin
        using .FiniteDiff
        include("finite_diff.jl")
    end

    @require FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000" begin
        using .FiniteDifferences
        include("finite_differences.jl")
    end

    @require ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210" begin
        using .ForwardDiff
        include("forward_diff.jl")
    end

    @require Manifolds = "1cead3c2-87b3-11e9-0ccd-23c62b72b94e" begin
        using .Manifolds
        include("manifolds.jl")
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
