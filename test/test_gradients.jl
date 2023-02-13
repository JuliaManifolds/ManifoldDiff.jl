using Manifolds
using ManifoldsBase

@testset "Gradients" begin
    M = Sphere(2)
    p = [0.0, 0.0, 1.0]
    q = [0.0, 1.0, 0.0]
    r = [1.0, 0.0, 0.0]
    @testset "Gradient of the distance function" begin
        X = zero_vector(M, p)
        ManifoldDiff.grad_distance!(M, X, p, q)
        Y = ManifoldDiff.grad_distance(M, p, q)
        Z = [0.0, 0.0, -π / 2] # known solution
        @test X == Y
        @test X == Z
        U = zero_vector(M, p)
        ManifoldDiff.grad_distance!(M, U, p, q, 1)
        V = ManifoldDiff.grad_distance(M, p, q, 1)
        W = -distance(M, q, p)^(-1) * log(M, q, p) # solution
        @test U == V
        @test U == W
        w = q
        U = zero_vector(M, w)
        ManifoldDiff.grad_distance!(M, U, w, q, 1)
        V = ManifoldDiff.grad_distance(M, w, q, 1)
        W = zero_vector(M, w) # solution
        @test U == V
        @test U == W
    end
end
