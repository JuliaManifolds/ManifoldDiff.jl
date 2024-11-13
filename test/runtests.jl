
using Test
using LinearAlgebra
using ManifoldDiff
using ForwardDiff
using FiniteDifferences
using Manifolds
using ManifoldsBase
using RecursiveArrayTools

@testset "ManifoldDiff.jl" begin
    include("test_differentials.jl")
    include("test_adjoint_differentials.jl")
    include("differentiation.jl")
    include("manifold_specializations.jl")
    include("test_gradients.jl")
    include("test_derivatives.jl")
    include("test_proximal_maps.jl")
    include("test_subgradients.jl")
    include("test_jacobians.jl")
end
