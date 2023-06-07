using Manifolds
using ManifoldsBase
using Random

@testset "Subgradients" begin
    M = Sphere(2)
    p = [0.0, 0.0, 1.0]
    @testset "Subgradient of the distance function" begin
        q = p
        X = zero_vector(M, q)
        Random.seed!(42)
        ManifoldDiff.subgrad_distance!(M, X, q, p)
        Y = ManifoldDiff.subgrad_distance(M, q, p)
        @test is_vector(M, p, X)
        @test norm(M, p, X) ≤ 1.0
        @test is_vector(M, p, Y)
        @test norm(M, p, Y) ≤ 1.0
    end
end
