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
    Weingarten,
    allocate_result,
    get_iterator,
    _write,
    _read

using Requires

using ADTypes:
    AbstractADType,
    AutoFiniteDiff,
    AutoFiniteDifferences,
    AutoForwardDiff,
    AutoReverseDiff,
    AutoZygote
import DifferentiationInterface as DI

struct NoneDiffBackend end

"""
    _derivative(f, t[, backend])

Compute the derivative of a callable `f` at time `t` computed using the given `backend`. If the backend is not explicitly
specified, it is obtained using the function [`default_differential_backend`](@ref).

This function calculates plain Euclidean derivatives, for Riemannian differentiation see
for example [`differential`](@ref differential(::AbstractManifold, ::Any, ::Real, ::AbstractRiemannianDiffBackend)).

!!! note

    Not specifying the backend explicitly will usually result in a type instability
    and decreased performance.
"""
function _derivative end

_derivative(f, t) = _derivative(f, t, default_differential_backend())

function _derivative!(f, X, t, backend = default_differential_backend())
    return copyto!(X, _derivative(f, t, backend))
end

"""
    _gradient(f, p[, backend])

Compute the gradient of a callable `f` at point `p` computed using the given `backend`. If the backend is not explicitly
specified, it is obtained using the function [`default_differential_backend`](@ref).

This function calculates plain Euclidean gradients, for Riemannian gradient calculation see
for example [`gradient`](@ref gradient(::AbstractManifold, ::Any, ::Any, ::AbstractRiemannianDiffBackend)).

!!! note

    Not specifying the backend explicitly will usually result in a type instability
    and decreased performance.
"""
function _gradient end

_gradient(f, p) = _gradient(f, p, default_differential_backend())

function _gradient!(f, X, p, backend = default_differential_backend())
    return copyto!(X, _gradient(f, p, backend))
end

"""
    _hessian(f, p[, backend])

Compute the Hessian of a callable `f` at point `p` computed using the given `backend`. If the backend is not explicitly
specified, it is obtained using the function [`default_differential_backend`](@ref).

This function calculates plain Euclidean Hessian.

!!! note

    Not specifying the backend explicitly will usually result in a type instability
    and decreased performance.
"""
function _hessian end

_hessian(f, p) = _hessian(f, p, default_differential_backend())

"""
    _jacobian(f, p[, backend])

Compute the Jacobian of a callable `f` at point `p` computed using the given `backend`. If the backend is not explicitly
specified, it is obtained using the function [`default_differential_backend`](@ref).

This function calculates plain Euclidean Jacobians, for Riemannian Jacobian calculation see
for example [`gradient`](@ref gradient(::AbstractManifold, ::Any, ::Any, ::AbstractRiemannianDiffBackend)).

!!! note

    Not specifying the backend explicitly will usually result in a type instability
    and decreased performance.
"""
function _jacobian end

_jacobian(f, p) = _jacobian(f, p, default_differential_backend())

function _jacobian!(f, X, p, backend = default_differential_backend())
    return copyto!(X, _jacobian(f, p, backend))
end

"""
    CurrentDiffBackend(backend)

A mutable struct for storing the current differentiation backend in a global
constant [`_current_default_differential_backend`](@ref).

# See also

[`default_differential_backend`](@ref), [`set_default_differential_backend!`](@ref)
"""
mutable struct CurrentDiffBackend
    backend::Any
end

"""
    _current_default_differential_backend

The instance of [`CurrentDiffBackend`](@ref) that stores the globally default
differentiation backend.
"""
const _current_default_differential_backend = CurrentDiffBackend(NoneDiffBackend())
"""
    default_differential_backend()

Get the default differentiation backend.
"""
default_differential_backend() = _current_default_differential_backend.backend

"""
    set_default_differential_backend!(backend)

Set current backend for differentiation to `backend`.
"""
function set_default_differential_backend!(backend)
    _current_default_differential_backend.backend = backend
    return backend
end

include("differentiationinterface.jl")

include("diagonalizing_projectors.jl")

include("adjoint_differentials.jl")
include("derivatives.jl")
include("differentials.jl")
include("gradients.jl")
include("Jacobi_fields.jl")
include("proximal_maps.jl")
include("subgradients.jl")
include("jacobians.jl")

include("riemannian_diff.jl")
include("embedded_diff.jl")



function __init__()
    # There is likely no way to set defaults without Requires.jl
    @require FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000" begin
        using .FiniteDifferences: central_fdm
        if default_differential_backend() === NoneDiffBackend()
            set_default_differential_backend!(AutoFiniteDifferences(central_fdm(5, 1)))
        end
    end

    @require ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210" begin
        if default_differential_backend() === NoneDiffBackend()
            set_default_differential_backend!(AutoForwardDiff())
        end
    end

    @require FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41" begin
        if default_differential_backend() === NoneDiffBackend()
            set_default_differential_backend!(AutoFiniteDiff())
        end
    end

    @require ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267" begin
        if default_differential_backend() === NoneDiffBackend()
            set_default_differential_backend!(AutoReverseDiff())
        end
    end

    @require Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f" begin
        if default_differential_backend() === NoneDiffBackend()
            set_default_differential_backend!(AutoZygote())
        end
    end

    return nothing
end

export riemannian_gradient,
    riemannian_gradient!,
    riemannian_Hessian,
    riemannian_Hessian!,
    subgrad_distance,
    subgrad_distance!
export set_default_differential_backend!, default_differential_backend

export TangentDiffBackend, RiemannianProjectionBackend

end # module
