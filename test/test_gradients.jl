using Manifolds
using ManifoldsBase

import ManifoldsBase: Weingarten!

Weingarten!(::Sphere, Y, p::Array, X, V) = (Y .= -X * p' * V)

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
    @testset "Riemannian Gradient / Hessian Conversion" begin
        M = Sphere(2)
        E = ℝ^3
        A = [1.0 0.1 0.2; 0.1 2.0 0.3; 0.2 0.3 3.0]
        ef(E, p) = 1 / 2 * p'A * p
        grad_ef(E, p) = A * p
        grad_f(M, p) = (A * p - (p'A * p) * p)
        Hess_ef(E, p, X) = A * X
        Hess_f(M, p, X) = A * X - (p' * A * X) .* p - (p' * A * p) .* X

        p = [1.0, 0.0, 0.0]
        X = [0.0, 0.1, 0.2]

        Y1 = riemannian_gradient(M, p, grad_ef(E, p))
        Y2 = grad_f(M, p)
        @test isapprox(M, p, Y1, Y2)

        Z1 = riemannian_Hessian(M, p, grad_ef(E, p), Hess_ef(M, p, X), X)
        Z2 = Hess_f(M, p, X)
        @test isapprox(M, p, Z1, Z2)
    end
end
