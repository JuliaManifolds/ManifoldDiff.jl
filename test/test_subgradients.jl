using Manifolds
using ManifoldsBase
using Random

@testset "Subgradients" begin
    M = Sphere(2)
    p = [0.0, 0.0, 1.0]
    q = [0.0, 1.0, 0.0]
    r = [1.0, 0.0, 0.0]
    @testset "Subgradient of the distance function" begin
        w = p
        X = zero_vector(M, w)
        Random.seed!(42)
        ManifoldDiff.subgrad_distance!(M, X, w, p)
        Y = ManifoldDiff.subgrad_distance(M, w, p)
        @test is_vector(M, p, X)
        @test norm(M, p, X) ≤ 1.0
        @test is_vector(M, p, Y)
        @test norm(M, p, Y) ≤ 1.0
        Z = ManifoldDiff.grad_distance(M, p, r)
        W = ManifoldDiff.subgrad_distance(M, p, r, 2)
        @test Z == W

        X = zero_vector(M, p)
        ManifoldDiff.subgrad_distance!(M, X, p, q, 2)
        Y = ManifoldDiff.subgrad_distance(M, p, q, 2)
        Z = [0.0, 0.0, -π / 2] # known solution
        @test X == Y
        @test X == Z
        U = zero_vector(M, p)
        ManifoldDiff.subgrad_distance!(M, U, p, q)
        V = ManifoldDiff.subgrad_distance(M, p, q)
        W = -distance(M, q, p)^(-1) * log(M, q, p) # solution
        @test U == V
        @test U == W
    end
end
